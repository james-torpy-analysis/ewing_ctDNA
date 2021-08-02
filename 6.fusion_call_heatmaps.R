
#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
ref_dir <- paste0(project_dir, "refs/")

func_dir <- paste0(project_dir, "scripts/functions/5.fusion_call_heatmaps/")
fusion_dir <- paste0(project_dir, "results/fusions/")
VAF_dir <- paste0(project_dir, "results/VAF_calculation/")

detection_dir <- paste0(project_dir, "results/detection_heatmaps/")
Robject_dir <- paste0(detection_dir, "/Rdata/")
plot_dir <- paste0(VAF_dir, "plots/")

system(paste0("mkdir -p ", func_dir))
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

library(dplyr)


####################################################################################
### 0. Load functions and colours ###
####################################################################################

fetch_fusion_no <- dget(paste0(func_dir, "fetch_fusion_no.R"))
longitudinal_patient_heatmap <- dget(
  paste0(func_dir, "longitudinal_patient_heatmap.R")
)
longitudinal_dilution_heatmap <- dget(
  paste0(func_dir, "longitudinal_dilution_heatmap.R")
)

hm_cols <- c(
  pathology_detection = "#75EA3D",
  no_pathology_detection = "#58B9DB",
  stringent_detection = "#F7A006",
  less_stringent_detection = "#F4D30B",
  supporting_reads = "white",
  no_supporting_reads = "black",
  unknown = "grey"
)

annot_cols <- c(
  FP = "#991425",
  VAF = "#221699",
  sread = "#115E0F"
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

# fetch VAFs with most supporting reads for each sample, 
# prioritising high conf:
VAFs <- unlist(
  lapply(samplenames, function(x) {
    
    print(x)
    
    all_VAFs <- tryCatch(
      all_VAFs <- readRDS(paste0(VAF_dir, x, "/Rdata/VAFs.Rdata")),
      error=function(err) NA
    )
    
    if (suppressWarnings(!is.na(all_VAFs))) {
      return(max(all_VAFs$VAF_fwd))
    } else {
      return(0)
    }
    
  })
)

VAF_df <- data.frame(
  VAF = VAFs
)
VAF_df$Sample <- gsub(
  "_.*$", "", 
  sub(".*?_(.+)", "\\1", unlist(samplenames))
)

# merge VAFs with fusion df:
fusion_dfs <- lapply(fusion_dfs, function(x) {
  x <- merge(
    x,
    VAF_df,
    by = "Sample"
  )
  return(x)
})

# save VAFs:
VAF_df <- VAF_df %>% column_to_rownames("Sample")
saveRDS(VAF_df, paste0(Robject_dir, "all_VAFs.Rdata"))

# fetch fusion-supporting reads:
fusion_read_nos <- sapply(samplenames, function(x) {

  print(x)

  return(
    read.table(
      paste0(VAF_dir, x, "/non_specific_fusion_supporting_reads.tsv"),
      header = T
    )[1,1]
  )
})
fusion_read_nos <- data.frame(
  Sample = gsub(
    "_.*$", "",
    sub(".*?_(.+)", "\\1", unlist(samplenames))
  ),
  Supporting_read_pairs = fusion_read_nos
)

# add to fusion_dfs:
fusion_dfs <- lapply(fusion_dfs, function(x) {
  return(
    merge(
      x,
      fusion_read_nos,
      by = "Sample"
    )
  )
})

save.image(paste0(Robject_dir, "temp_img.Rdata"))
#load(paste0(Robject_dir, "temp_img.Rdata"))


####################################################################################
### 2. Create longitudinal heatmaps ###
####################################################################################

patient_heatmaps_VAF <- longitudinal_patient_heatmap(
  fusion_df = fusion_dfs$patient,
  hm_title = "Patient EWSR1/FLI1 fusion detections",
  annot = "VAF",
  hm_cols = hm_cols
)

dilution_heatmaps_VAF <- longitudinal_dilution_heatmap(
  fusion_df = fusion_dfs$dilution,
  hm_title = "Cell line EWSR1/FLI1 fusion detections",
  annot = "VAF",
  hm_cols = hm_cols
)

all_heatmaps <- list(
  patient_heatmap_FP = longitudinal_patient_heatmap(
    fusion_df = fusion_dfs$patient,
    hm_title = "Patient EWSR1/FLI1 fusion detections",
    annot = "false positives",
    hm_cols = hm_cols
  ),
  patient_heatmap_VAF = patient_heatmaps_VAF$VAF_annot,
  patient_heatmap_sread = patient_heatmaps_VAF$sread_annot,
  dilution_heatmap_FP = longitudinal_dilution_heatmap(
    fusion_df = fusion_dfs$dilution,
    hm_title = "Cell line EWSR1/FLI1 fusion detections",
    annot = "false positives",
    hm_cols = hm_cols
  ),
  dilution_heatmap_VAF = dilution_heatmaps_VAF$VAF_annot,
  dilution_heatmap_sread = dilution_heatmaps_VAF$sread_annot
)

for (i in seq_along(all_heatmaps)) {
  
  # convert to grob:
  hm_grob <- grid.grabExpr(
    draw(all_heatmaps[[i]], gap = unit(6, "mm"))
  )
  dev.off()
  
  if (length(grep("patient", names(all_heatmaps)[i])) > 0) {
    
    png(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.png"),
      width = 10,
      height = 6,
      unit = "in",
      res = 300
    )
    annot_coord <-  0.41
    
  } else {
    
    png(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.png"),
      width = 10,
      height = 3,
      unit = "in",
      res = 300
    )
    annot_coord <-  0.37
    
  }
  
    grid.newpage()
    
    pushViewport(viewport(x = 0.027, y = 0.01, width = 0.964, height = 0.89, 
                          just = c("left", "bottom")))
      grid.draw(hm_grob)
    popViewport()
    
    if (length(grep("FP", names(all_heatmaps)[i])) > 0) {
      
      pushViewport(viewport(x = 0.99, y = annot_coord, width = 0.2, height = 0.1, 
        just = c("right", "top")))
      
        grid.text(
          "x = No. false positives", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "FP"])
        )
      
      popViewport()
      
    } else if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      
      pushViewport(viewport(x = 0.95, y = annot_coord, width = 0.2, height = 0.1, 
        just = c("right", "top")))
        grid.text(
          "x% = VAF", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "VAF"])
        )
      popViewport()
      
    } else {
      
      pushViewport(viewport(x = 0.995, y = annot_coord, width = 0.2, height = 0.1, 
        just = c("right", "top")))
        grid.text(
          "x = No. supporting reads", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "sread"])
        )
      popViewport()
    }
  
  dev.off()
  
}

