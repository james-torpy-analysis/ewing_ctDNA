
args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
#projectname <- "ewing_ctDNA"
#samplename <- "409_012_combined" 

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
result_dir <- paste0(project_dir, "results/")
bam_dir <- paste0(result_dir, "picard/")
fusion_dir <- paste0(result_dir, "fusions/", samplename, "/")
out_path <- paste0(result_dir, "VAF_calculation/", samplename, "/less_stringent/")

Robject_dir <- paste0(out_path, "Rdata/")
hist_dir <- paste0(out_path, "histograms/")
out_bam_dir <- paste0(out_path, "bams/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", hist_dir))
system(paste0("mkdir -p ", out_bam_dir))

func_dir <- paste0(project_dir, "scripts/functions/")

library(GenomicAlignments)
library(Rsamtools)
library(scales)

source(paste0(func_dir, "vaf_functions.R"))

min_overlap_R1 <- 19
min_overlap_R2 <- 19

# define fusion grs and orientations:
fusion_gr <- list(
  EWSR1_FLI1=list(
    EWSR1=GRanges(
      seqnames="chr22",
      ranges=IRanges(start=29664257, end=29696511),
      strand="*" ),
    FLI1=GRanges(
      seqnames="chr11",
      ranges=IRanges(start=128556430, end=128683162),
      strand="*" )),
  EWSR1_ETV1=list(
    EWSR1=GRanges(
      seqnames="chr22",
      ranges=IRanges(start=29664257, end=29696511),
      strand="*" ),
    ETV1=GRanges(
      seqnames="chr7",
      ranges=IRanges(start=13930854, end=14031050),
      strand="*" )),
  EWSR1_ERG=list(
    EWSR1=GRanges(
      seqnames="chr22",
      ranges=IRanges(start=29664257, end=29696511),
      strand="*" ),
    ERG=GRanges(
      seqnames="chr21",
      ranges=IRanges(start=39739183, end=40033707),
      strand="*" )))

fusion_orientations <- c(EWSR1_FLI1="+", EWSR1_ETV1="+", EWSR1_ERG="+")


## 1) read bam file

file_bam <- paste0(
    bam_dir,
    samplename, "/", samplename, ".dedup.sorted.by.coord.bam" )

what <- c(
    "qname",
    "flag",
    "rname",
    "strand",
    "pos",
    "qwidth",
    "mapq",
    "cigar",
    "mrnm",
    "mpos",
    "isize",
    "seq",
    "qual" )

param <- ScanBamParam(
    flag = scanBamFlag(isUnmappedQuery = F),
    what = what)

ga <- readGAlignments(file_bam, use.names = TRUE, param = param)
gr <- granges(ga, use.mcols = TRUE)

## check all reads paired
stopifnot(all(mcols(gr)$flag %% 2 >= 1))

# ## remove reads with segments not properly aligned
# remove <- unique(names(gr)[which(mcols(gr)$flag %% 4 < 2)])
# gr <- gr[which(!names(gr) %in% remove), ]

## check no unmapped reads
stopifnot(all(mcols(gr)$flag %% 8 < 4))

## check no multi-mappers
stopifnot(all(mcols(gr)$flag %% 512 < 256))

## add whether alignment is for first or second read in template
mcols(gr)$R1 <- mcols(gr)$flag %% 128 >= 64
mcols(gr)$R2 <- mcols(gr)$flag %% 256 >= 128
stopifnot(all(mcols(gr)$R1 + mcols(gr)$R2 == 1L))

## extract R1 + R2
## R1 ~ gene-specific primer
## R2 ~ universal primer
tmp  <- gr[mcols(gr)$R1]
mcols(tmp) <- NULL
R1 <- split(tmp, names(tmp))
tmp  <- gr[mcols(gr)$R2]
mcols(tmp) <- NULL
R2 <- split(tmp, names(tmp))

# save as RDS:
saveRDS(gr, paste0(Robject_dir, "filtered_reads.Rdata"))


## 2) detect fusions

for (i in seq_along(fusion_gr)) {

  print(paste0("Detecting ", names(fusion_gr)[i], " fusions..."))

  # find fusion supporting split reads:
  split_R1 <- R1[lengths(range(R1)) == 2L]
  split_R2 <- R2[lengths(range(R2)) == 2L]
  
#  # require that primary and supp alignments are on the same strand, to filter out
#  # weird fusions e.g. upstream of gene 1 joined to upstream of gene 2:
#  split_R1 <- split_R1[sapply(split_R1, function(x) length(unique(strand(x))) == 1)]
#  split_R2 <- split_R2[sapply(split_R2, function(x) length(unique(strand(x))) == 1)]

  # find fusion supporting split R1s:
  fusion_R1 <- sapply(split_R1, function (x) {
    i1 <- which(x %over% fusion_gr[[i]][[1]])
    i2 <- which(x %over% fusion_gr[[i]][[2]])
    if (length(i1) == 0 || length(i2) == 0) {
      return()
    } else {
      paste(as.character(flank(x[i1], 1, FALSE)), 
        as.character(flank(x[i2], 1, TRUE)) )}})
  fusion_R1 <- unlist(fusion_R1)
  
  # find fusion supporting split R2s:
  fusion_R2 <- sapply(split_R2, function (x) {
    i1 <- which(x %over% fusion_gr[[i]][[1]])
    i2 <- which(x %over% fusion_gr[[i]][[2]])
    if (length(i1) == 0 || length(i1) == 0) {
      return()
    } else {
      paste(as.character(flank(x[i1], 1, TRUE)),
        as.character(flank(x[i2], 1, FALSE)) )}})
  fusion_R2 <- unlist(fusion_R2)
  
  # save as bams:
  try(
    writeSam(file_bam, names(fusion_R1), paste0(out_bam_dir, "/", names(fusion_gr)[i], 
      "_fusion_split_R1s.sam" )))
  try(
    writeSam(file_bam, names(fusion_R2), paste0(out_bam_dir, "/", names(fusion_gr)[i], 
      "_fusion_split_R2s.sam" )))

  # keep only unique fusions:
  fusion_list <- as.list(unique(c(gsub(":\\+|:\\-", ":\\*", fusion_R1), 
    gsub(":\\+|:\\-", ":\\*", fusion_R2)) ))
  fusions <- lapply(fusion_list, function(x) {
    spl <- strsplit(x, " ")[[1]]
    fusion_gr <- as(spl[1], "GRanges")
    spl2 <- strsplit(spl[2], ":")[[1]]
    fusion_gr$join_chr <- spl2[1]
    fusion_gr$join_coord <- as.numeric(spl2[2])
    return(fusion_gr) })
  fusions <- unlist(as(fusions, "GRangesList"))
  
  if (length(fusions) >= 1) {
    
    ## add orientation for identified translocations
    strand(fusions) <- fusion_orientations[i]
    mcols(fusions)$join_strand <- fusion_orientations[i]
    
    for (j in seq_along(fusions)) {
      
      fusion <- fusions[j]
  
      gene_a_breakpoint <- fusion    
      gene_b_breakpoint <- GRanges(
        fusion$join_chr,
        IRanges(
          fusion$join_coord,
          fusion$join_coord),
        fusion$join_strand)
      mcols(gene_a_breakpoint) <- NULL
      
      breakpoint <- gsub(":", "_", as.character(gene_a_breakpoint), fixed = TRUE)
      print(breakpoint)
      
      ## consider 1Mb upstream or downstream of gene 1 breakpoint
      gene_a_upstream <- flank(gene_a_breakpoint, 1e6, start = TRUE)
      gene_a_upstream <- resize(gene_a_upstream, 1e6 + 1) ## add breakpoint position
      gene_a_dnstream <- flank(gene_a_breakpoint, 1e6, start = FALSE)
      gene_a_dnstream_rev <- reverseStrand(gene_a_dnstream)
      
      ## consider 1Mb upstream or downstream of fusion partner breakpoint
      gene_b_upstream <- flank(gene_b_breakpoint, 1e6, start = TRUE)
      gene_b_upstream <- resize(gene_b_upstream, 1e6 + 1) ## add breakpoint position
      gene_b_dnstream <- flank(gene_b_breakpoint, 1e6, start = FALSE)
      gene_b_dnstream_rev <- reverseStrand(gene_b_dnstream)
      
      ## non-supporting reads that satisfy overlap criteria for driver fusion
      nonsupp_fwd <- intersect(
        names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_a_dnstream_rev))) >= min_overlap_R2)))
      ## supporting reads that satisfy overlap criteria for driver fusion
      supp_fwd <- intersect(
        names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_b_dnstream_rev))) >= min_overlap_R2)))
      
      ## diagnost plots for overlap
      R1_nonsupp_fwd <- R1[which(names(R1) %in% nonsupp_fwd)]
      R2_nonsupp_fwd <- R2[which(names(R2) %in% nonsupp_fwd)]
      R1_supp_fwd <- R1[which(names(R1) %in% supp_fwd)]
      R2_supp_fwd <- R2[which(names(R2) %in% supp_fwd)]
      
      pdf(paste0(hist_dir, "hist_overlap_", breakpoint, "_fwd.pdf"))
      par(mfrow = c(2, 2))
      if (length(pintersect(R1_nonsupp_fwd, gene_a_upstream))>0) {
        hist(-sum(width(pintersect(R1_nonsupp_fwd, gene_a_upstream))),
          xlim = c(-180, 0), xlab = "Overlap [bp]",
          main = paste0("fusion non-supporting ", names(fusion_gr[[i]])[1], " upstream") )}
      if (length(pintersect(R2_nonsupp_fwd, gene_a_dnstream_rev))>0) {
        hist(sum(width(pintersect(R2_nonsupp_fwd, gene_a_dnstream_rev))),
           xlim = c(0, 180), xlab = "Overlap [bp]",
           main = paste0("fusion non-supporting ", names(fusion_gr[[i]])[1], " dnstream") )}
      if (length(pintersect(R1_supp_fwd, gene_a_upstream))>0) {
        hist(-sum(width(pintersect(R1_supp_fwd, gene_a_upstream))),
             xlim = c(-180, 0), xlab = "Overlap [bp]",
             main = paste0("fusion supporting ", names(fusion_gr[[i]])[1], " upstream") )}
      if (length(pintersect(R1_supp_fwd, gene_a_upstream))>0) {
        hist(sum(width(pintersect(R1_supp_fwd, gene_a_upstream))),
             xlim = c(0, 180), xlab = "Overlap [bp]",
             main = paste0("fusion supporting ", names(fusion_gr[[i]])[2], " dnstream") )}
      dev.off()
      
      ## non-supporting reads that satisfy overlap criteria for reciprocal fusion
      nonsupp_rev <- intersect(
        names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_a_upstream))) >= min_overlap_R2)))
      ## supporting reads that satisfy overlap criteria for reciprocal fusion
      supp_rev <- intersect(
        names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_b_upstream))) >= min_overlap_R2)))
      
      ## diagnost plots for overlap
      R1_nonsupp_rev <- R1[which(names(R1) %in% nonsupp_rev)]
      R2_nonsupp_rev <- R2[which(names(R2) %in% nonsupp_rev)]
      R1_supp_rev <- R1[which(names(R1) %in% supp_rev)]
      R2_supp_rev <- R2[which(names(R2) %in% supp_rev)]
      
      pdf(paste0(hist_dir, "hist_overlap_", breakpoint, "_rev.pdf"))
      par(mfrow = c(2, 2))
      if (length(pintersect(R2_nonsupp_rev, gene_a_upstream))>0) {
        hist(-sum(width(pintersect(R2_nonsupp_rev, gene_a_upstream))),
          xlim = c(-180, 0), xlab = "Overlap [bp]",
          main = paste0("reciproc non-supporting ", names(fusion_gr[[i]])[1], " upstream") )}
      if (length(pintersect(R1_nonsupp_rev, gene_a_dnstream_rev))>0) {
        hist(sum(width(pintersect(R1_nonsupp_rev, gene_a_dnstream_rev))),
          xlim = c(0, 180), xlab = "Overlap [bp]",
          main = paste0("reciproc non-supporting ", names(fusion_gr[[i]])[1], " dnstream") )}
      if (length(pintersect(R2_supp_rev, gene_b_upstream))>0) {
        hist(-sum(width(pintersect(R2_supp_rev, gene_b_upstream))),
            xlim = c(-180, 0), xlab = "Overlap [bp]",
            main = paste0("reciproc supporting ", names(fusion_gr[[i]])[2], " upstream"), 
            silent = TRUE )}
      if (length(pintersect(R1_supp_rev, gene_a_dnstream_rev))>0) {
        hist(sum(width(pintersect(R1_supp_rev, gene_a_dnstream_rev))),
            xlim = c(0, 180), xlab = "Overlap [bp]",
            main = paste0("reciproc supporting ", names(fusion_gr[[i]])[1], " dnstream") , 
            silent = TRUE )}
      dev.off()
      
      ## calculate VAFs for both translocations
      print(length(nonsupp_fwd))
      print(length(supp_fwd))
      VAF_fwd <- length(supp_fwd) / (length(nonsupp_fwd) + length(supp_fwd))
      print(VAF_fwd)
      
      print(length(nonsupp_rev))
      print(length(supp_rev))
      VAF_rev <- length(supp_rev) / (length(nonsupp_rev) + length(supp_rev))
      print(VAF_rev)
      
      if (j==1) {
        VAFs <- list(
          data.frame(VAF_fwd, VAF_rev)
        )
      } else {
        VAFs[[j]] <- data.frame(VAF_fwd, VAF_rev)
      }
      names(VAFs)[j] <- paste0("fusion_", start(fusions[j]))
  
      VAFs[[j]]$Fusion <- paste0(
        gsub(":\\+", "-", as.character(fusion)), fusion$join_chr, ":", fusion$join_coord )
      VAFs[[j]]$Forward_supporting <- length(supp_fwd)
      VAFs[[j]]$Forward_non_supporting <- length(nonsupp_fwd)
      VAFs[[j]]$Forward_total <- VAFs[[j]]$Forward_supporting + VAFs[[j]]$Forward_non_supporting
      VAFs[[j]]$Reverse_supporting <- length(supp_rev)
      VAFs[[j]]$Reverse_non_supporting <- length(nonsupp_rev)
      VAFs[[j]]$Reverse_total <- VAFs[[j]]$Reverse_supporting + VAFs[[j]]$Reverse_non_supporting
      
      ## write reads to SAM for inspection
      writeSam(file_bam, nonsupp_fwd, paste0(out_bam_dir, "/reads_", breakpoint, "_nonsupp_fwd.sam"))
      writeSam(file_bam, supp_fwd, paste0(out_bam_dir, "/reads_", breakpoint, "_supp_fwd.sam"))
      writeSam(file_bam, nonsupp_rev, paste0(out_bam_dir, "/reads_", breakpoint, "_nonsupp_rev.sam"))
      writeSam(file_bam, supp_rev, paste0(out_bam_dir, "/reads_", breakpoint, "_supp_rev.sam"))
      
    }
    
    # collate VAFs as df:
    VAF_df <- do.call("rbind", VAFs)
  
    # remove reverse fusions for now:
    VAF_df <- subset(VAF_df, 
      select = -c(VAF_rev, Reverse_supporting, Reverse_non_supporting, Reverse_total) )
  
    # remove fusions with no supporting reads:
    VAF_df <- VAF_df[VAF_df$Forward_supporting != 0,]
  
    if (nrow(VAF_df) > 0) {
      # order by forward supporting read number:
      VAF_df <- VAF_df[order(VAF_df$Forward_supporting),]

      # add fusion type column and bind dfs:
      VAF_df$Fusion_type <- names(fusion_gr)[i]
      if (exists("all_VAF")) {
        all_VAF <- rbind(all_VAF, VAF_df)
      } else {
        all_VAF <- VAF_df 
      }
    }
  }
}

if (exists("all_VAF")) {
  write.table(
    all_VAF,
    paste0(out_path, "VAF.tsv"),
    sep = "\t",
    quote = F,
    row.names = T,
    col.names = T)
  saveRDS(all_VAF, file.path(Robject_dir, "VAF.rds"))
} else {
  ## create output file for snakemake
  all_VAF <- NULL
  saveRDS(all_VAF, file.path(Robject_dir, "VAF.rds"))
}
