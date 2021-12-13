#!/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/envs/ewing_ctDNA/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
projectname <- "ewing_ctDNA"
samplename <- "409_004_combined"

home_dir <- "/share/ScratchGeneral/jamtor"
#home_dir <- "/Users/torpor/clusterHome"
project_dir <- file.path(home_dir, "projects", projectname)
func_dir <- file.path(project_dir, "scripts/functions")
result_dir <- file.path(project_dir, "results")
bam_dir <- file.path(result_dir, "picard")
fusion_dir <- file.path(result_dir, "fusions", samplename)
out_path <- file.path(result_dir, "VAF_calculation", samplename)

Robject_dir <- file.path(out_path, "Rdata")
hist_dir <- file.path(out_path, "histograms")
out_bam_dir <- file.path(out_path, "bams")

dir.create(Robject_dir, recursive=T)
dir.create(hist_dir, recursive=T)
dir.create(out_bam_dir, recursive=T)

func_dir <- file.path(project_dir, "scripts/functions/")

libs <- c("GenomicAlignments", "Rsamtools", "scales")
lapply(libs, library, character.only=T)

source(file.path(func_dir, "vaf_functions.R"))

min_overlap_R1 <- 19
min_overlap_R2 <- 19
max_suppl_dist <- 10  # supplementary reads must fall within this no. nt for read
# to be counted as supporting
min_suppl_reads <- 2  # minimum supplementary (split) reads required to call 
# a supporting read in absence of discordant (non-split) alignments

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
gene_orients <- c(EWSR1_FLI1="++", EWSR1_ETV1="+-", EWSR1_ERG="+-")

## 1) read bam file

file_bam <- file.path(bam_dir,samplename, 
  paste0(samplename, ".dedup.sorted.by.coord.bam") )

print(paste0("Reading and filtering bam file ", file_bam, "..."))

what <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar",
  "mrnm", "mpos", "isize", "seq","qual" )

param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F), what = what)

ga <- readGAlignments(file_bam, use.names = TRUE, param = param)
gr <- granges(ga, use.mcols = TRUE)

## check all reads paired
stopifnot(all(mcols(gr)$flag %% 2 >= 1))

## check no unmapped reads
stopifnot(all(mcols(gr)$flag %% 8 < 4))

## remove reads with unmapped mates
gr <- gr[mcols(gr)$flag %% 16 < 8]
stopifnot(all(mcols(gr)$flag %% 16 < 8))

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
mcols(tmp) <- subset(mcols(tmp), select=flag)
R1 <- split(tmp, names(tmp))
tmp  <- gr[mcols(gr)$R2]
mcols(tmp) <- subset(mcols(tmp), select=flag)
R2 <- split(tmp, names(tmp))

# save as RDS:
saveRDS(gr, file.path(Robject_dir, "filtered_reads.Rdata"))


## 2) detect fusions

print("Detecting all fusions...")

# find fusion supporting split reads:
split_R1 <- R1[lengths(range(R1)) == 2L]
split_R2 <- R2[lengths(range(R2)) == 2L]
saveRDS(list(split_R1=split_R1, split_R2=split_R2), 
  file.path(Robject_dir, "split_reads.rds"))

# save as bams:
writeSam(file_bam, names(split_R1), file.path(out_bam_dir, "split_R1s.sam"))
writeSam(file_bam, names(split_R2), file.path(out_bam_dir, "split_R2s.sam"))

for (j in seq_along(fusion_gr)) {

  print(paste0("Detecting ", names(fusion_gr)[j], " fusions..."))

  # set up cluster for parLapply:
  cl <- makePSOCKcluster(7)
  clusterExport(cl, c("fusion_gr", "j"))
  clusterEvalQ(cl, library(GenomicAlignments))

  # find fusion supporting split R1s:
  print("Finding fusion supporting R1s...")
  fusion_R1 <- parSapply(cl, split_R1, function (x) {
    i1 <- which(x %over% fusion_gr[[j]][[1]])
    i2 <- which(x %over% fusion_gr[[j]][[2]])
    if (length(i1) == 0 || length(i2) == 0) {
      return()
    } else {
      paste(as.character(flank(x[i1], 1, FALSE)), as.character(flank(x[i2], 1, TRUE)) )}})
  fusion_R1 <- unlist(fusion_R1)
  
  # find fusion supporting split R2s:
  print("Finding fusion supporting R2s...")
  fusion_R2 <- parSapply(cl, split_R2, function (x) {
    i1 <- which(x %over% fusion_gr[[j]][[1]])
    i2 <- which(x %over% fusion_gr[[j]][[2]])
    if (length(i1) == 0 || length(i2) == 0) {
      return()
    } else {
      paste(as.character(flank(x[i1], 1, TRUE)), as.character(flank(x[i2], 1, FALSE)) )}})
  fusion_R2 <- unlist(fusion_R2)

  stopCluster(cl)
    
  # keep only unique fusions:
  print("Removing fusion entry duplicates...")
  fusion_list <- as.list(unique(c(gsub(":\\+|:\\-", ":\\*", fusion_R1), gsub(":\\+|:\\-", ":\\*", fusion_R2))))
  fusions <- lapply(fusion_list, function(x) {
    spl <- strsplit(x, " ")[[1]]
    fusion_gr <- as(spl[1], "GRanges")
    spl2 <- strsplit(spl[2], ":")[[1]]
    fusion_gr$join_chr <- spl2[1]
    fusion_gr$join_coord <- as.numeric(spl2[2])
    return(fusion_gr) })
  fusions <- unlist(as(fusions, "GRangesList"))

  if (j==1) {
    all_fusions <- list(fusions)
  } else {
    all_fusions[[j]] <- fusions }
}
names(all_fusions) <- names(fusion_gr)


## 3) filter fusions and calculate VAFs

print("Filtering fusions and calculating VAFs...")

# filter fusions and calculate VAFs:
for (j in seq_along(all_fusions)) {
  if (length(all_fusions[[j]]) >= 1) {
    ## add orientation for identified translocations
    strand(all_fusions[[j]]) <- substring(gene_orients[j], 1, 1)
    mcols(all_fusions[[j]])$join_strand <- substring(gene_orients[j], 2, 2)

    # filter fusions and calculate VAFs:
    fgenes <- strsplit(names(all_fusions)[j], "_")[[1]]
    VAFs <- filter_and_VAF(fusions_gr=all_fusions[[j]], max_suppl_dist, min_suppl_reads,
      gene_a_name=fgenes[1], gene_b_name=fgenes[2], R1, R2, hist_dir, out_bam_dir )

    # collate VAFs as df:
    VAF_df <- do.call("rbind", VAFs)
    # remove fusions with no supporting reads:
    VAF_df <- VAF_df[VAF_df$All_supporting != 0,]
    if (nrow(VAF_df) > 0) {
      # order by forward supporting read number:
      VAF_df <- VAF_df[order(VAF_df$All_supporting),]
      # add fusion type column and bind dfs:
      VAF_df$Fusion_type <- names(all_fusions)[j]
      if (exists("all_VAF")) {
        all_VAF <- rbind(all_VAF, VAF_df)
        # reorder VAF df:
        all_VAF <- subset(all_VAF, select=c("Fusion_type", "Fusion", "VAF_forward", 
          "VAF_up_to_upstream", "VAF_reverse", "VAF_down_to_downstream", 
          "Forward_non_supporting", "Forward_supporting", "Forward_total", 
          "Up_to_upstream_supporting", "Up_to_upstream_total", "Reverse_non_supporting", 
          "Down_to_upstream_supporting", "Down_to_upstream_total", 
          "Down_to_downstream_supporting", "Down_to_downstream_total", "All_supporting" ))
      } else {all_VAF <- VAF_df } }}}

print("Saving VAFs...")

if (exists("all_VAF")) {
  write.table(all_VAF, file.path(out_path, "VAF.tsv"), sep="\t", quote=F, 
    row.names=T, col.names=T)
  saveRDS(all_VAF, file.path(Robject_dir, "VAF.rds"))
} else { saveRDS(NULL, file.path(Robject_dir, "VAF.rds")) } # for snakemake
