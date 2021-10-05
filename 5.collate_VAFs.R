projectname <- "ewing_ctDNA"

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome"
project_dir <- file.path(home_dir, "projects", projectname)
func_dir <- file.path(project_dir, "scripts/functions")
result_dir <- file.path(project_dir, "results/")
ref_dir <- file.path(project_dir, "refs")
variant_dir <- file.path(result_dir, "smcounter2/tables")
in_path <- file.path(result_dir, "VAF_calculation")

Robject_dir <- file.path(in_path, "Rdata")
table_dir <- file.path(in_path, "tables")
plot_dir <- file.path(in_path, "plots")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", plot_dir))

library(ggplot2)
library(cowplot)
library(ggrepel)

compare_VAF <- dget(file.path(func_dir, "compare_VAF.R"))


####################################################################################
### 1. Load VAFs ###
####################################################################################

summary_df <- read.table(
  file.path(variant_dir, "sample_summary_with_smcounter_SNV.tsv"),
  sep = "\t",
  header = TRUE )

VAFs <- lapply(summary_df$Library_id, function(x) {
  print(x)
  VAF <- as.data.frame(readRDS(file.path(in_path, x, "Rdata/VAF.rds")))
  if (!is.null(VAF)) {
    if (nrow(VAF) > 0) {
      VAF <- VAF[order(VAF$Forward_supporting, decreasing = T),]
      VAF$Library_id <- rep(x, nrow(VAF))
      VAF$VAF_fwd <- round(as.numeric(VAF$VAF_fwd), 4)
    } else {
      VAF <- data.frame(
        VAF_fwd = "not_detected",
        Fusion = "not_detected",
        Forward_supporting = "not_detected",
        Forward_non_supporting = "not_detected",
        Forward_total = "not_detected",
        Library_id = x
      )
      rownames(VAF) <- x
    }
  }
  return(VAF)
})

VAF_df <- do.call("rbind", VAFs)
VAF_df <- subset(VAF_df, select = c(Library_id, Fusion, VAF_fwd,
  Forward_supporting, Forward_non_supporting, Forward_total ) )


####################################################################################
### 2. Add previous VAF estimates and write table ###
####################################################################################

# format and make values numeric:
all_VAFs <- merge(summary_df, VAF_df, by="Library_id")

all_VAFs <- subset(all_VAFs, select = c(
  Patient_id, Sample_id, Library_id, Site, Treatment.dilution, 
  Sanger_TP53_SNV, TP53_SNV_type, smCounter2_TP53_SNV, ddPCR_TP53_VAF, 
  GeneGlobe_TP53_VAF, smCounter2_TP53_VAF, smCounter2_TP53_effect_size,
  smCounter2_TP53_UMT, smCounter2_TP53_VMT, smCounter2_TP53_qual, 
  Sanger_STAG2_SNV, smCounter2_STAG2_SNV, ddPCR_STAG2_VAF, 
  GeneGlobe_STAG2_VAF, smCounter2_STAG2_VAF, smCounter2_STAG2_effect_size,
  smCounter2_STAG2_UMT, smCounter2_STAG2_VMT, smCounter2_STAG2_qual, 
  Pathology_EWSR1_FLI1, Fusion, VAF_fwd,
  Forward_supporting, Forward_non_supporting, Forward_total ))

# convert fusion VAFs to percentages:
all_VAFs$VAF_fwd[all_VAFs$VAF_fwd != "not_detected"] <- 
  as.numeric(all_VAFs$VAF_fwd[all_VAFs$VAF_fwd != "not_detected"])*100

colnames(all_VAFs) <- c(
  "Patient_id", "Sample_id", "Library_id", "Site", "Treatment/dilution", 
  "Sanger_TP53_SNV", "TP53_SNV_type", "smCounter2_TP53_SNV", "ddPCR_TP53_VAF", 
  "GeneGlobe_TP53_VAF", "smCounter2_TP53_VAF", "smCounter2_TP53_effect_size",
  "smCounter2_TP53_UMT", "smCounter2_TP53_VMT", "smCounter2_TP53_qual", 
  "Sanger_STAG2_SNV", "smCounter2_STAG2_SNV", "ddPCR_STAG2_VAF", 
  "GeneGlobe_STAG2_VAF", "smCounter2_STAG2_VAF", "smCounter2_STAG2_effect_size", 
  "smCounter2_STAG2_UMT", "smCounter2_STAG2_VMT", "smCounter2_STAG2_qual", 
  "Pathology_EWSR1_FLI1", "Fusion_EWSR1_FLI1", "Fusion_VAF",
  "Forward_supporting", "Forward_non_supporting", "Forward_total" )

# subset deletion columns only, order by patient and write:
fusion_VAFs <- subset(all_VAFs, select = -c(
  Sanger_TP53_SNV, TP53_SNV_type, smCounter2_TP53_SNV, 
  ddPCR_TP53_VAF, GeneGlobe_TP53_VAF, smCounter2_TP53_VAF, 
  smCounter2_TP53_effect_size,
  smCounter2_TP53_UMT, smCounter2_TP53_VMT, smCounter2_TP53_qual, 
  Sanger_STAG2_SNV, smCounter2_STAG2_SNV, ddPCR_STAG2_VAF, 
  GeneGlobe_STAG2_VAF, smCounter2_STAG2_VAF, smCounter2_STAG2_effect_size,
  smCounter2_STAG2_UMT, smCounter2_STAG2_VMT, smCounter2_STAG2_qual ))
fusion_VAFs$Library_id <- factor(
  fusion_VAFs$Library_id, levels = summary_df$Library_id )
fusion_VAFs <- fusion_VAFs[order(fusion_VAFs$Library_id),]

write.table(
  fusion_VAFs, 
  file.path(table_dir, "final_fusion_summary.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE )

# remove secondary deletion rows, order by patient and write:
all_VAFs <- all_VAFs[!duplicated(all_VAFs$Library_id),]
all_VAFs <- all_VAFs[match(summary_df$Library_id, all_VAFs$Library_id), ]
write.table(
  all_VAFs, 
  file.path(table_dir, "final_summary.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE )


####################################################################################
### 3. Plot VAF correlations ###
####################################################################################

# check distributions of VAFs:
plot(hist(as.numeric(all_VAFs$smCounter2_TP53_VAF[
  all_VAFs$smCounter2_TP53_VAF != 0 & all_VAFs$smCounter2_TP53_VAF != "not_detected" ])))
plot(hist(as.numeric(all_VAFs$smCounter2_STAG2_VAF[
  all_VAFs$smCounter2_STAG2_VAF != 0 & all_VAFs$smCounter2_STAG2_VAF != "not_detected" ])))
dev.off()

# create smCounter2 TP53 vs STAG2 VAF plot:
VAF_df <- subset(all_VAFs, select = c(
  Sample_id, Treatment.dilution, smCounter2_TP53_VAF, smCounter2_STAG2_VAF ))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment[43:54] <- "ES8"
VAF_df$treatment[55:65] <- "A673"
TP53_vs_STAG2_VAF <- compare_VAF(
  VAF_df, lab1 = "smCounter2_TP53", lab2 = "smCounter2_STAG2", 
  lim = 100, cortype = "pearson" )

pdf(file.path(plot_dir, "smCounter2_TP53_vs_STAG2_VAFs.pdf"),
    height = 3.5,
    width = 5)
  print(TP53_vs_STAG2_VAF)
dev.off()

png(file.path(plot_dir, "smCounter2_TP53_vs_STAG2_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(TP53_vs_STAG2_VAF)
dev.off()

# create smCounter2 TP53 vs fusion VAF plot:
VAF_df <- subset(all_VAFs, select = c(
  Sample_id, Treatment, smCounter2_TP53_VAF, Fusion_VAF ))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment[43:54] <- "ES8"
VAF_df$treatment[55:65] <- "A673"

TP53_vs_fusion_VAF <- compare_VAF(
  VAF_df, lab1 = "smCounter2_TP53", lab2 = "Fusion", 
  lim = 100, cortype = "pearson" )

pdf(file.path(plot_dir, "smCounter2_TP53_vs_fusion_VAFs.pdf"),
    height = 3.5,
    width = 5)
print(TP53_vs_fusion_VAF)
dev.off()

png(file.path(plot_dir, "smCounter2_TP53_vs_fusion_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(TP53_vs_fusion_VAF)
dev.off()

# create smCounter2 STAG2 vs fusion VAF plot:
VAF_df <- subset(all_VAFs, select = c(
  Sample_id, Treatment.dilution, smCounter2_STAG2_VAF, Fusion_VAF ))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment[43:54] <- "ES8"
VAF_df$treatment[55:65] <- "A673"

STAG2_vs_fusion_VAF <- compare_VAF(
  VAF_df, lab1 = "smCounter2_STAG2", lab2 = "Fusion", 
  lim = 100, cortype = "pearson" )

pdf(file.path(plot_dir, "smCounter2_STAG2_vs_fusion_VAFs.pdf"),
    height = 3.5,
    width = 5)
print(STAG2_vs_fusion_VAF)
dev.off()

png(file.path(plot_dir, "smCounter2_STAG2_vs_fusion_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(STAG2_vs_fusion_VAF)
dev.off()

