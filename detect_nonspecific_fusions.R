
args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
projectname <- "ewing_ctDNA"
samplename <- "409_009_combined" 

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

# define fusion grs and orientations:
EWSR1_gr=GRanges(
  seqnames="chr22",
  ranges=IRanges(start=29664257, end=29696511),
  strand="*" )

fusion_orientations <- c(EWSR1="+")


################################################################################
### 1. Load split reads ###
################################################################################

split_reads <- readRDS(file.path(Robject_dir, "split_reads.rds"))
split_R1 <- split_reads[[1]]
split_R2 <- split_reads[[2]]




for (j in seq_along(fusion_gr)) {
  # set up cluster for parLapply:
  cl <- makePSOCKcluster(7)
  clusterExport(cl, c("fusion_gr", "j"))
  clusterEvalQ(cl, library(GenomicAlignments))

  # find fusion supporting split R1s:
  fusion_R1 <- parSapply(cl, split_R1, function (x) {
    i1 <- which(x %over% fusion_gr[[j]][[1]])
    i2 <- which(x %over% fusion_gr[[j]][[2]])
    if (length(i1) == 0 || length(i2) == 0) {
        return()
    } else {
      paste(as.character(flank(x[i1], 1, FALSE)), 
        as.character(flank(x[i2], 1, TRUE)) )}})
  fusion_R1 <- unlist(fusion_R1)
  
  # find fusion supporting split R2s:
  fusion_R2 <- parSapply(cl, split_R2, function (x) {
    i1 <- which(x %over% fusion_gr[[j]][[1]])
    i2 <- which(x %over% fusion_gr[[j]][[2]])
    if (length(i1) == 0 || length(i2) == 0) {
        return()
    } else {
      paste(as.character(flank(x[i1], 1, TRUE)),
          as.character(flank(x[i2], 1, FALSE)) )}})
  fusion_R2 <- unlist(fusion_R2)

  stopCluster(cl)
    
  # keep only unique fusions:
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


## 2) filter fusions and calculate VAFs

print("Filtering fusions and calculating VAFs...")

# filter fusions and calculate VAFs:
for (j in seq_along(all_fusions)) {
  if (length(all_fusions[[j]]) >= 1) {
    ## add orientation for identified translocations
    strand(all_fusions[[j]]) <- fusion_orientations[1]
    mcols(all_fusions[[j]])$join_strand <- fusion_orientations[1]

    # filter fusions and calculate VAFs:
    fgenes <- strsplit(names(all_fusions)[j], "_")[[1]]
    VAFs <- filter_and_VAF(fusions_gr=all_fusions[[j]], gene_a_name = fgenes[1], 
      gene_b_name = fgenes[2], hist_dir, out_bam_dir )
    
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
      VAF_df$Fusion_type <- names(all_fusions)[j]
      if (exists("all_VAF")) {
        all_VAF <- rbind(all_VAF, VAF_df)
      } else {
        all_VAF <- VAF_df }}}}

print("Saving VAFs...")

if (exists("all_VAF")) {
  write.table(all_VAF, paste0(out_path, "VAF.tsv"), sep = "\t", quote = F, 
    row.names = T, col.names = T )
  saveRDS(all_VAF, file.path(Robject_dir, "VAF.rds"))
} else {
  ## create output file for snakemake
  all_VAF <- NULL
  saveRDS(all_VAF, file.path(Robject_dir, "VAF.rds")) }
