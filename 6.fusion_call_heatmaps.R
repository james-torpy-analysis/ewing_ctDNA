
home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
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
library(naturalsort)


####################################################################################
### 0. Load functions and colours ###
####################################################################################

fetch_fusion_no <- dget(paste0(func_dir, "fetch_fusion_no.R"))
mutation_heatmap <- dget(paste0(func_dir, "mutation_heatmap.R") )

hm_cols <- c(
  pathology_detection = "#75EA3D",
  no_pathology_detection = "#58B9DB",
  detected = "#F4D30B",
  not_detected = "black",
  no_sample = "grey" )

annot_cols <- c(
  VAF = "#991425",
  sread = "#430F82" )


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
  levels = c("tumour", "naive", "NACT1", "NACT2", "resection", "ACT1", 
    "targeted", "relapse", "ACT2" ))
patient_df <- patient_df[order(patient_df$Treatment.dilution),]

# remove STAG2 detection in ES_6 ACT2:
patient_df[
  patient_df$Patient_id == "ES_6" & 
  patient_df$Treatment.dilution == "ACT2",]$smCounter2_STAG2_VAF <- "not_detected"

# separate and order cell line samples:
cl_df <- summary_df[summary_df$Site == "cell_line" & !is.na(summary_df$Site),]
cl_df <- cl_df[naturalorder(cl_df$Treatment.dilution, decreasing = TRUE),]
cl_df <- cl_df[order(cl_df$Patient_id),]


####################################################################################
### 2. Create longitudinal heatmaps ###
####################################################################################

patient_heatmaps <- mutation_heatmap(
  fusion_df = patient_df,
  type = "patient",
  hm_cols )

cl_heatmaps <- mutation_heatmap(
  fusion_df = cl_df,
  type = "dilution",
  hm_cols )

all_heatmaps <- c(patient_heatmaps, cl_heatmaps)
names(all_heatmaps) <- c(
  paste0("patient_heatmap_", names(all_heatmaps)[1:2]),
  paste0("cell_line_heatmap_", names(all_heatmaps)[3:4]) )

for (i in seq_along(all_heatmaps)) {
  
  # convert to grob:
  hm_grob <- grid.grabExpr(
    draw(all_heatmaps[[i]], gap = unit(6, "mm"))
  )
  dev.off()
  
  # write to png:
  if (length(grep("patient", names(all_heatmaps)[i])) > 0) {
    
    png(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.png"),
      width = 14.5,
      height = 9,
      unit = "in",
      res = 300
    )
    annot_y <-  0.43
    
    if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      annot_x <- 0.798
    } else {
      annot_x <- 0.828
    }
    
  } else {
    
    png(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.png"),
      width = 20,
      height = 3,
      unit = "in",
      res = 300
    )
    annot_y <-  0.41
    
    if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      annot_x <- 0.834
    } else {
      annot_x <- 0.83
    }
    
  }
  
    grid.newpage()
    
    pushViewport(viewport(x = 0.027, y = 0.01, width = 0.964, height = 0.89, 
                          just = c("left", "bottom")))
      grid.draw(hm_grob)
    popViewport()
    
    if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      
      pushViewport(viewport(x = annot_x, y = annot_y, width = 0.2, height = 0.1, 
        just = c("left", "top")))
        grid.text(
          "x% = VAF", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "VAF"])
        )
      popViewport()
      
    } else {
      
      pushViewport(viewport(x = annot_x, y = annot_y, width = 0.2, height = 0.1, 
        just = c("left", "top")))
        grid.text(
          "x = No. supporting reads", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "sread"])
        )
      popViewport()
    }
  
  dev.off()
  
  # write to pdf:
  if (length(grep("patient", names(all_heatmaps)[i])) > 0) {
    
    pdf(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.pdf"),
      width = 14.5,
      height = 9,
    )
    annot_y <-  0.41
    
  } else {
    
    pdf(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.pdf"),
      width = 26,
      height = 3,
    )
    annot_y <-  0.37
    
  }
  
  grid.newpage()
  
  pushViewport(viewport(x = 0.027, y = 0.01, width = 0.964, height = 0.89, 
                        just = c("left", "bottom")))
  grid.draw(hm_grob)
  popViewport()
  
  if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
    
    pushViewport(viewport(x = annot_x, y = annot_y, width = 0.2, height = 0.1, 
                          just = c("left", "top")))
    grid.text(
      "x% = VAF", 
      gp=gpar(fontsize=10, fontface="bold", 
              col=annot_cols[names(annot_cols) == "VAF"])
    )
    popViewport()
    
  } else {
    
    pushViewport(viewport(x = annot_x, y = annot_y, width = 0.2, height = 0.1, 
                          just = c("left", "top")))
    grid.text(
      "x = No. supporting reads", 
      gp=gpar(fontsize=10, fontface="bold", 
              col=annot_cols[names(annot_cols) == "sread"])
    )
    popViewport()
  }
  
  dev.off()
  
}

