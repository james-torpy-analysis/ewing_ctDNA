
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

dir.create(func_dir)
dir.create(Robject_dir)
dir.create(plot_dir)

library(dplyr)
library(tibble)


####################################################################################
### 0. Load functions and colours ###
####################################################################################

fetch_fusion_no <- dget(paste0(func_dir, "fetch_fusion_no.R"))
longitudinal_patient_heatmap <- dget(
  paste0(func_dir, "longitudinal_patient_heatmap.R") )
longitudinal_dilution_heatmap <- dget(
  paste0(func_dir, "longitudinal_dilution_heatmap.R") )

hm_cols <- c(
  pathology_detection = "#75EA3D",
  no_pathology_detection = "#58B9DB",
  detected = "#F7A006",
  not_detected = "black",
  unknown = "grey" )

annot_cols <- c(
  FP = "#991425",
  VAF = "#221699",
  sread = "#115E0F" )


####################################################################################
### 1. Load summary table and split into patient and cell line dfs ###
####################################################################################

summary_df <- read.table(
  file.path(VAF_dir, "tables/final_summary.tsv"),
  header = T,
  sep = "\t" )

# separate and order patients:
patient_df <- summary_df[
  summary_df$Site != "cell_line" | is.na(summary_df$Site), ]

patient_df$Treatment.dilution <- factor(patient_df$Treatment.dilution,
  levels = c("tumour", "resection", "naive", "NACT1", "NACT2", "ACT1", 
    "targeted", "relapse", "ACT2" ))
patient_df <- patient_df[order(patient_df$Treatment.dilution),]

# separate and order cell line samples:
cl_df <- summary_df[summary_df$Site == "cell_line",]
cl_df <- cl_df[order(cl_df$Treatment.dilution, decreasing = TRUE),]


####################################################################################
### 2. Create longitudinal heatmaps ###
####################################################################################

patient_heatmaps_VAF <- longitudinal_patient_heatmap(
  fusion_df = patient_df,
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

