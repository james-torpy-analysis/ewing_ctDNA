hm_cols <- c(
  FISH_yes = "#75EA3D",
  FISH_no = "#D68EB7",
  yes = "#F4D30B", 
  no = "black", 
  unknown = "grey"
)
hm_na_col <- "grey"

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
in_path <- paste0(project_dir, "results/svaba/BWA_and_picard/")
non_collapsed_path <- paste0(project_dir, "results/svaba/non_collapsed/")

func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/")
ref_dir <- paste0(project_dir, "refs/")
Robject_dir <- paste0(project_dir, "results/fusions/")

system(paste0("mkdir -p ", Robject_dir))


####################################################################################
### 0. Load packages, fusions and GOI coordinates ###
####################################################################################

library(rtracklayer)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(knitr)

load_breakpoints <- dget(paste0(func_dir, "load_breakpoints.R"))
find_EWSR1_FLI1_fusions <- dget(paste0(func_dir, "find_EWSR1_FLI1_fusions.R"))
create_summary_tables <- dget(paste0(func_dir, "create_summary_tables.R"))
longitudinal_heatmap <- dget(
  paste0(func_dir, "longitudinal_heatmap.R")
)

if (!file.exists(paste0(Robject_dir, "EWSR1_GOI_fusions.Rdata"))) {
  # define GOI co-ordinates:
  GOI <- list(
    EWSR1 = GRanges(
      seqnames = "chr22",
      ranges = IRanges(start = 29663998, end = 29696515),
      strand = "*"
    ),
    FLI1 = GRanges(
      seqnames = "chr11",
      ranges = IRanges(start = 128554430, end = 128685162),
      strand = "*"
    ),
    ERG = GRanges(
      seqnames = "chr7",
      ranges = IRanges(start = 13928854, end = 14033050),
      strand = "*"
    ),
    ETV1 = GRanges(
      seqnames = "chr21",
      ranges = IRanges(start = 39737183, end = 40035704),
      strand = "*"
    ),
    ETV4 = GRanges(
      seqnames = "chr17",
      ranges = IRanges(start = 41603214, end = 41625708),
      strand = "*"
    ),
    FEV = GRanges(
      seqnames = "chr2",
      ranges = IRanges(start = 219843809, end = 219851906),
      strand = "*"
    )
  )
  
  # load metadata:
  patient_df <- read.table(
    paste0(ref_dir, "ES_samples_by_patient.tsv"),
    sep = "\t",
    header = T
  )
  # format colnames and Sample column:
  colnames(patient_df) <- gsub("\\.", "_", colnames(patient_df))
  samplenames <- patient_df$Sample
  patient_df$Sample <- gsub(
    "_.*$", "", 
    gsub("409_", "", patient_df$Sample)
  )
  
  # order rows:
  split_df <- split(patient_df, patient_df$Patient)
  
  split_df <- lapply(split_df, function(x) {
    
    # tumour samples first:
    df <- x[x$Treatment == "tumour",]
    # add resected:
    df <- rbind(df, x[x$Treatment == "resection",])
    # add treatment naive:
    df <- rbind(df,x[x$Treatment == "naive",])
    # add NACT:
    df <- rbind(df, x[grep("NACT", x$Treatment),])
    # add pre-relapse ACT:
    df <- rbind(df,x[x$Treatment == "ACT1",])
    # add targeted:
    df <- rbind(df,x[x$Treatment == "targeted",])
    # add relapses:
    df <- rbind(df, x[x$Treatment == "relapse",])
    # add post-relapse ACT:
    return(rbind(df,x[x$Treatment == "ACT2",]))
    
  })
  # merge split elements:
  patient_df <- do.call("rbind", split_df)
  
  # load dilution samples metadata:
  dilution_df <- read.table(
    paste0(ref_dir, "ES_samples_by_dilution.tsv"),
    sep = "\t",
    header = T
  )
  # order rows and format Sample column:
  colnames(dilution_df) <- gsub("\\.", "_", colnames(dilution_df))
  dilution_df <- arrange(dilution_df, desc(as.numeric(Dilution)))
  samplenames <- c(samplenames, dilution_df$Sample)
  samplenames <- sort(samplenames)
  dilution_df$Sample <- gsub(
    "_.*$", "", 
    gsub("409_", "", dilution_df$Sample)
  )
  
  
  ####################################################################################
  ### 1. Load VCFs and define EWSR1 co-ordinates ###
  ####################################################################################
  
  bp <- list(
    collapsed = list(
      high_conf_bp = load_breakpoints(samplenames, in_path),
      low_conf_bp = load_breakpoints(samplenames, in_path, high_conf = FALSE)
    ),
    non_collapsed = list(
      high_conf_bp = load_breakpoints(samplenames, non_collapsed_path),
      low_conf_bp = load_breakpoints(samplenames, non_collapsed_path, high_conf = FALSE)
    )
  )
  
  
  ####################################################################################
  ### 2. Find overlaps with breakpoints and EWSR1 ###
  ####################################################################################
  
  fusion_results <- lapply(bp, function(x) {
    lapply(
      x,
      find_EWSR1_FLI1_fusions,
      GOI_list = GOI,
      patient_df, 
      dilution_df
    )
  })
    
  
  # save EWSR1 fusions:
  saveRDS(
    fusion_results,
    paste0(Robject_dir, "EWSR1_GOI_fusions.Rdata")
  )
  
} else {
  fusion_results <- readRDS(paste0(Robject_dir, "EWSR1_GOI_fusions.Rdata"))
}


####################################################################################
### 3. Create summary tables ###
####################################################################################

summary_tables <- lapply(fusion_results, function(x) {
  create_summary_tables(x)
})

####################################################################################
### 3. Create longitudinal heatmaps ###
####################################################################################

# patient heatmaps:
patient_detections <- lapply(fusion_results, function(x) {
  lapply(x, function(y) y$fusion_nos$patient_df)
})

for (i in 1:length(patient_detections)) {
  if (i==1) {
    patient_heatmaps <- list(
      lapply(
        patient_detections[[i]], 
        longitudinal_heatmap,
        hm_title = paste0(
          "UMI ", 
          gsub("_", " ", names(patient_detections)[i])
        ),
        type = "patient",
        hm_cols = hm_cols
      )
    )
  } else {
    patient_heatmaps[[i]] <- lapply(
      patient_detections[[i]], 
      longitudinal_heatmap,
      hm_title = paste0(
        "UMI ", 
        gsub("_", "-", names(patient_detections)[i])
      ),
      type = "patient",
      hm_cols = hm_cols
    )
  }
}
names(patient_heatmaps) <- names(patient_detections)

# dilution heatmaps:
dilution_detections <- lapply(fusion_results, function(x) {
  lapply(x, function(y) y$fusion_nos$dilution_df)
})

for (i in 1:length(dilution_detections)) {
  if (i==1) {
    dilution_heatmaps <- list(
      lapply(
        dilution_detections[[i]], 
        longitudinal_heatmap,
        hm_title = paste0(
          "UMI ", 
          gsub("_", " ", names(dilution_detections)[i])
        ),
        type = "dilution",
        hm_cols = hm_cols
      )
    )
  } else {
    dilution_heatmaps[[i]] <- lapply(
      dilution_detections[[i]], 
      longitudinal_heatmap,
      hm_title = paste0(
        "UMI ", 
        gsub("_", "-", names(dilution_detections)[i])
      ),
      type = "dilution",
      hm_cols = hm_cols
    )
  }
}
names(dilution_heatmaps) <- names(dilution_detections)

