
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

func_dir <- paste0(project_dir, "scripts/functions/")
fusion_dir <- paste0(project_dir, "results/fusions/")
VAF_dir <- paste0(project_dir, "results/VAF_calculation/")


####################################################################################
### 0. Load patient data, fusion calls and VAFs ###
####################################################################################

######
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
######
samplenames <- list.files(VAF_dir, pattern = "409")

fusion_results <- readRDS()

VAFs <- sapply(samplenames, function(x) {
  
})





####################################################################################
### 3. Create summary tables ###
####################################################################################

summary_tables <- lapply(fusion_results, function(x) {
  create_summary_tables(x)
})


####################################################################################
### 1. Create longitudinal heatmaps ###
####################################################################################

# patient heatmaps with false positive annotation:
patient_detections <- lapply(fusion_results, function(x) {
  lapply(x, function(y) y$fusion_nos$patient_df)
})

for (i in 1:length(patient_detections)) {
  if (i==1) {
    patient_heatmaps <- list(
      list(
        false_annot = lapply(
          patient_detections[[i]], 
          longitudinal_heatmap,
          hm_title = paste0(
            "UMI ", 
            gsub("_", " ", names(patient_detections)[i])
          ),
          type = "patient",
          hm_cols = hm_cols,
          SNP_annot = FALSE
        ),
        short_variant_annot = lapply(
          patient_detections[[i]], 
          longitudinal_heatmap,
          hm_title = paste0(
            "UMI ", 
            gsub("_", " ", names(patient_detections)[i])
          ),
          type = "patient",
          hm_cols = hm_cols,
          SNP_annot = TRUE
        )
      )
    )
  } else {
    patient_heatmaps[[i]] <- list(
      false_annot = lapply(
        patient_detections[[i]], 
        longitudinal_heatmap,
        hm_title = paste0(
          "UMI ", 
          gsub("_", " ", names(patient_detections)[i])
        ),
        type = "patient",
        hm_cols = hm_cols,
        SNP_annot = FALSE
      ),
      short_variant_annot = lapply(
        patient_detections[[i]], 
        longitudinal_heatmap,
        hm_title = paste0(
          "UMI ", 
          gsub("_", " ", names(patient_detections)[i])
        ),
        type = "patient",
        hm_cols = hm_cols,
        SNP_annot = TRUE
      )
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
