
args = commandArgs(trailingOnly=TRUE)

samplename <- args[1]
#samplename <- "409_027_DCKVC_TAGGCATG-CTCTCTAT_L001"

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
in_dir <- paste0(project_dir, "results/svaba/BWA_and_picard/", samplename, "/")
non_collapsed_dir <- paste0(project_dir, "results/svaba/non_collapsed/", samplename, "/")

func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/")
ref_dir <- paste0(project_dir, "refs/")
Robject_dir <- paste0(project_dir, "results/fusions/", samplename, "/")

system(paste0("mkdir -p ", Robject_dir))


####################################################################################
### 0. Load packages, fusions and GOI coordinates ###
####################################################################################

library(rtracklayer)
library(dplyr)
library(tibble)
library(ComplexHeatmap)

load_breakpoints <- dget(paste0(func_dir, "load_breakpoints.R"))
find_EWSR1_FLI1_fusions <- dget(paste0(func_dir, "find_EWSR1_FLI1_fusions.R"))
longitudinal_heatmap <- dget(
  paste0(func_dir, "longitudinal_heatmap.R")
)

#if (!file.exists(paste0(Robject_dir, "EWSR1_GOI_fusions.Rdata"))) {
  
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
  
  
  ####################################################################################
  ### 1. Load VCFs and define EWSR1 co-ordinates ###
  ####################################################################################
  
  bp <- list(
    high_conf_bp = load_breakpoints(samplename, in_dir),
    low_conf_bp = load_breakpoints(samplename, in_dir, high_conf = FALSE)
  )
  
  
  ####################################################################################
  ### 2. Find overlaps with breakpoints and EWSR1 ###
  ####################################################################################
  
  fusion_results <- lapply(
    bp,
    find_EWSR1_FLI1_fusions,
    GOI_list = GOI
  )
  
  # save EWSR1 fusions:
  saveRDS(
    fusion_results,
    paste0(Robject_dir, "EWSR1_GOI_fusions.Rdata")
  )
  
#}

