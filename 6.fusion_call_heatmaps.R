
#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
ref_dir <- paste0(project_dir, "refs/")

func_dir <- paste0(project_dir, "scripts/functions/")
fusion_dir <- paste0(project_dir, "results/fusions/")
VAF_dir <- paste0(project_dir, "results/VAF_calculation/")

library(dplyr)


####################################################################################
### 0. Load functions and colours ###
####################################################################################

fetch_fusion_no <- dget(paste0(func_dir, "fetch_fusion_no.R"))

hm_cols <- c(
  FISH_yes = "#75EA3D",
  FISH_no = "#D68EB7",
  stringent_yes = "#F7A006",
  less_stringent_yes = "#F4D30B", 
  no = "black",
  unknown = "grey"
)


####################################################################################
### 1. Load patient data, fusion calls and VAFs ###
####################################################################################

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

# define sample names:
samplenames <- as.list(list.files(VAF_dir, pattern = "409"))

# fetch fusion numbers:
fusion_nos <- fetch_fusion_no(
  samplenames, 
  fusion_dir,
  patient_df,
  dilution_df
)

# add fusion numbers to metadata df:
for (i in 1:length(fusion_nos)) {
  
  if (i==1) {
    fusion_dfs <- list(
      patient = merge(patient_df, fusion_nos[[i]], by = "Sample")
    )
  } else {
    fusion_dfs$dilution = merge(dilution_df, fusion_nos[[i]], by = "Sample")
  }
  
}

VAFs <- unlist(
  lapply(samplenames, function(x) {
    tryCatch(
      read.table(
        paste0(VAF_dir, x, "/VAF.txt"),
        header = F
      ),
      error=function(err) NA
    )
  })
)

# turn NAs to 0:
VAFs[is.na(VAFs)] <- 0

VAF_df <- data.frame(
  VAF = VAFs
)
VAF_df$Sample <- gsub(
  "_.*$", "", 
  sub(".*?_(.+)", "\\1", unlist(samplenames))
)

fusion_dfs <- lapply(fusion_dfs, function(x) {
  x <- merge(
    x,
    VAF_df,
    by = "Sample"
  )
  return(x)
})

######

# fetch fusion-supporting reads:
fusion_read_nos <- sapply(samplenames, function(x) {
  
  tryCatch(
    read_no <- read.table(
      paste0(VAF_dir, x, "/final_fusion_read_nos.txt"),
      header = T,
      stringsAsFactors = F
    ),
    error=function(err) NA
  )
  
})

######


####################################################################################
### 2. Create longitudinal heatmaps ###
####################################################################################

patient_heatmap_FP <- longitudinal_heatmap(
  fusion_dfs$patient,
  hm_title = "Patient EWSR1/FLI1 fusion detections",
  type = "patient",
  annotation = "false positives",
  hm_cols = hm_cols
)

patient_heatmap_VAF <- longitudinal_heatmap(
  fusion_dfs$patient,
  hm_title = "Patient EWSR1/FLI1 fusion detections",
  type = "patient",
  annotation = "VAF",
  hm_cols = hm_cols
)

dilution_heatmap_FP <- longitudinal_heatmap(
  fusion_dfs$dilution,
  hm_title = "Cell line EWSR1/FLI1 fusion detections",
  type = "dilution",
  annotation = "false positives",
  hm_cols = hm_cols
)

dilution_heatmap_VAF <- longitudinal_heatmap(
  fusion_dfs$dilution,
  hm_title = "Cell line EWSR1/FLI1 fusion detections",
  type = "dilution",
  annotation = "VAF",
  hm_cols = hm_cols
)

