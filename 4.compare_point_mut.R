
home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")

variant_path <- paste0(project_dir, "results/smcounter2/")

plot_dir <- paste0(variant_path, "plots/")
table_dir <- paste0(variant_path, "tables/")
Robject_dir <- paste0(variant_path, "Rdata/")

system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", Robject_dir))

library(ggplot2)
library(rtracklayer)
library(GenomicAlignments)
library(tibble)
library(cowplot)
library(ggrepel)

fetch_sm_vafs <- dget(paste0(func_dir, "fetch_sm_vafs.R"))

meta <- read.table(paste0(ref_dir, "/metadata.tsv"), header = T)

# fetch sample and treatment info:
meta_list <- split(meta, meta$Library_id)

# define ROI:
roi_bed <- import(file.path(ref_dir, "CDHS-34925Z-409.roi.bed"))
TP53 <- roi_bed[seqnames(roi_bed) == "chr17"]
STAG2 <- roi_bed[seqnames(roi_bed) == "chrX"]

# define relevant germline variants:
gm_df <- data.frame(
  seqnames=c("chrX", "chrX"),
  start=c(123184949, 123195593),
  width=rep(1, 1),
  ref=c("C", "C"),
  alt=c("CT", "CT") )

gm_var <- GRanges(
  seqnames=gm_df$seqnames,
  ranges=IRanges(start=gm_df$start, width=gm_df$width),
  ref=gm_df$ref, alt=gm_df$alt )

# fetch and annotate VAFs:
sm_vars <- list(
  TP53 = do.call("rbind", lapply(meta_list, fetch_sm_vafs, roi=TP53, gm_var,
    reverse_strand=TRUE )),
  STAG2 = do.call("rbind", lapply(meta_list, fetch_sm_vafs, roi=STAG2, gm_var)) )

# format:
for (i in seq_along(sm_vars)) {
  sm_vars[[i]]$new_point_mut <- paste0(sm_vars[[i]]$seqnames, ":", sm_vars[[i]]$start, "_", 
    sm_vars[[i]]$ref, ">", sm_vars[[i]]$alt )
  var_df <- subset(sm_vars[[i]], select = c(new_point_mut, effect_size, VAF, UMT, VMT, qual) )
  colnames(var_df) <- c("smCounter2_point_mut", "smCounter2_effect_size", "smCounter2_VAF", 
    "smCounter2_UMT", "smCounter2_VMT", "smCounter2_qual" )
  colnames(var_df) <- gsub("2_", paste0("2_", names(sm_vars)[i], "_"), colnames(var_df))
  var_df$Library_id <- rownames(var_df)
  if (i==1) {
    var_dfs <- list(var_df)
  } else {
    var_dfs[[i]] <- var_df
  }
}

# merge with prior info:
final_vars <- merge(meta, var_dfs[[1]], by = "Library_id")
final_vars <- merge(final_vars, var_dfs[[2]], by = "Library_id")

# remove point_mut info for not_detected samples:
final_vars$smCounter2_TP53_point_mut[grep("not_detected", final_vars$smCounter2_TP53_point_mut)] <- 
  "not_detected"
final_vars$smCounter2_STAG2_point_mut[grep("not_detected", final_vars$smCounter2_STAG2_point_mut)] <- 
  "not_detected"

# subset:
final_df <- subset(final_vars, 
  select = c(Patient_id, Sample_id, Library_id, Site, Treatment.dilution, 
    Sanger_TP53_point_mut, smCounter2_TP53_point_mut, ddPCR_TP53_VAF, 
    GeneGlobe_TP53_VAF, smCounter2_TP53_VAF, smCounter2_TP53_effect_size, 
    smCounter2_TP53_UMT, smCounter2_TP53_VMT, smCounter2_TP53_qual, 
    Sanger_STAG2_point_mut, smCounter2_STAG2_point_mut, ddPCR_STAG2_VAF, 
    GeneGlobe_STAG2_VAF, smCounter2_STAG2_VAF, smCounter2_STAG2_effect_size, 
    smCounter2_STAG2_UMT, smCounter2_STAG2_VMT, smCounter2_STAG2_qual, 
    Pathology_EWSR1_FLI1 ) )

# sort by original metadata order:
final_df <- final_df[match(meta$Sample_id, final_df$Sample_id),]

write.table(
  final_df,
  file.path(table_dir, "sample_summary_with_smcounter_point_mut.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE )

