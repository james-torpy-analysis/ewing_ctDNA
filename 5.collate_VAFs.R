projectname <- "ewing_ctDNA"
germline_mutations <- c("chr17:7578457_G>A")

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome"
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

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggrepel)

compare_VAF <- dget(file.path(func_dir, "compare_VAF.R"))


####################################################################################
### 1. Load VAFs ###
####################################################################################

summary_df <- read.table(
  file.path(variant_dir, "sample_summary_with_smcounter_point_mut.tsv"),
  sep = "\t",
  header = TRUE )
summary_df <- subset(summary_df, select = c(
  Patient_id, Sample_id, Library_id, Site, Treatment.dilution, 
  Sanger_TP53_point_mut, smCounter2_TP53_point_mut, ddPCR_TP53_VAF, 
  GeneGlobe_TP53_VAF, smCounter2_TP53_VAF, smCounter2_TP53_effect_size,
  smCounter2_TP53_UMT, smCounter2_TP53_VMT, smCounter2_TP53_qual, 
  Sanger_STAG2_point_mut, smCounter2_STAG2_point_mut, ddPCR_STAG2_VAF, 
  GeneGlobe_STAG2_VAF, smCounter2_STAG2_VAF, smCounter2_STAG2_effect_size,
  smCounter2_STAG2_UMT, smCounter2_STAG2_VMT, smCounter2_STAG2_qual, 
  Pathology_EWSR1_FLI1, Pathology_EWSR1_ETV1, Pathology_EWSR1_ERG ))

VAFs <- lapply(summary_df$Library_id, function(x) {
  print(x)
  VAF <- as.data.frame(readRDS(file.path(in_path, x, "Rdata/VAF.rds")))
  if (!is.null(VAF)) {
    if (nrow(VAF) > 0) {
      VAF <- VAF[order(VAF$All_supporting, decreasing = T),]
      VAF$Library_id <- rep(x, nrow(VAF))
    } else {
      VAF <- data.frame(
        Fusion_type="not_detected",
        Fusion="not_detected",
        VAF_forward="not_detected",
        VAF_up_to_upstream="not_detected",
        VAF_reverse="not_detected",
        VAF_down_to_downstream="not_detected",
        Forward_non_supporting="not_detected",
        Forward_supporting="not_detected",
        Forward_total="not_detected",
        Up_to_upstream_supporting="not_detected",
        Up_to_upstream_total="not_detected",
        Reverse_non_supporting="not_detected",
        Down_to_upstream_supporting="not_detected",
        Down_to_upstream_total="not_detected",
        Down_to_downstream_supporting="not_detected",
        Down_to_downstream_total="not_detected",
        All_supporting="not_detected",
        Library_id = x
      )
      rownames(VAF) <- x
    }
  }
  return(VAF)
})
VAF_df <- do.call("rbind", VAFs)

# subset for only forward orientated fusions:
VAF_df <- subset(VAF_df, select = c(Library_id, Fusion_type, Fusion, VAF_forward,
  Forward_supporting, Forward_non_supporting, Forward_total ) )

# remove entries with no supporting reads:
VAF_df[VAF_df$Forward_supporting==0, grepl("Fusion|orward", colnames(VAF_df))] <- "not_detected"

# convert fusion VAFs to percentages:
VAF_df$VAF_forward[is.na(VAF_df$VAF_forward)] <- "not_detected"
VAF_df$VAF_forward[VAF_df$VAF_forward != "not_detected"] <- round(as.numeric(
  VAF_df$VAF_forward[VAF_df$VAF_forward != "not_detected"] )*100, 2)

# save for correlation scatterplots:
scatter_VAF <- VAF_df

# separate different fusions into different columns:
VAF_df <- VAF_df %>%
  mutate(id = row_number()) %>% 
  pivot_wider(
    names_from = Fusion_type,
    values_from = c(Fusion, VAF_forward, Forward_supporting, Forward_non_supporting, 
    Forward_total),
    names_glue = "{Fusion_type}_{.value}"
  ) %>% 
  select(-id) %>%
  select(-contains("not_detected")) %>%
  as.data.frame()
colnames(VAF_df) <- gsub("Fusion", "fusion", colnames(VAF_df))
colnames(VAF_df) <- gsub("Forward", "forward", colnames(VAF_df))

# change NAs to "not detected":
VAF_df[is.na(VAF_df)] <- "not_detected"

# order columns:
VAF_df <- VAF_df[,order(
  sapply(strsplit(colnames(VAF_df), "_"), function(x) x[2]), decreasing=T )]


####################################################################################
### 2. Add previous VAF estimates and write table ###
####################################################################################

# format and make values numeric:
all_VAFs <- merge(summary_df, VAF_df, by="Library_id")

# reorder to group fusion columns together:
for (i in grep("pathology", colnames(all_VAFs))) {
  fus <- gsub("_pathology", "", colnames(all_VAFs)[i])
  all_VAFs <- relocate(all_VAFs, colnames(all_VAFs)[i], 
    .before=colnames(all_VAFs)[grepl(fus, colnames(all_VAFs)) & 
    grepl("fusion", colnames(all_VAFs))]) }

# remove germline mutations:
all_VAFs[all_VAFs$smCounter2_TP53_point_mut %in% germline_mutations, 
  grepl("TP53", colnames(all_VAFs))] <- "not_detected"

# order rows by no. supporting reads, forward VAF:
all_VAFs$Library_id <- factor(all_VAFs$Library_id, levels=summary_df$Library_id)
all_VAFs <- arrange(all_VAFs, Library_id, desc(EWSR1_FLI1_forward_supporting), 
  desc(EWSR1_FLI1_VAF_forward) )
all_VAFs <- relocate(all_VAFs, c(Reads, UMIs, Reads_per_UMI), .after=last_col())

# subset fusion columns only, order by patient and write:
fusion_VAFs <- subset(all_VAFs, select = -c(
  Sanger_TP53_point_mut, smCounter2_TP53_point_mut, 
  ddPCR_TP53_VAF, GeneGlobe_TP53_VAF, smCounter2_TP53_VAF, 
  smCounter2_TP53_effect_size,
  smCounter2_TP53_UMT, smCounter2_TP53_VMT, smCounter2_TP53_qual, 
  Sanger_STAG2_point_mut, smCounter2_STAG2_point_mut, ddPCR_STAG2_VAF, 
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

# remove secondary deletion rows:
all_VAFs <- all_VAFs[!duplicated(all_VAFs$Library_id),]

# order by patient and write:
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

scatter_VAF <- merge(summary_df, scatter_VAF, by="Library_id")

colnames(scatter_VAF)[24] <- "EWSR1_FLI1_pathology"
colnames(scatter_VAF)[25] <- "EWSR1_ETV1_pathology"
colnames(scatter_VAF)[26] <- "EWSR1_ERG_pathology"

# order rows by no. supporting reads, forward VAF, then forward total:
#scatter_VAF <- arrange(scatter_VAF, desc(Forward_supporting), desc(VAF_forward))
# remove secondary deletion rows:
scatter_VAF <- scatter_VAF[!duplicated(scatter_VAF$Library_id),]

# check distributions of VAFs:
plot(hist(as.numeric(scatter_VAF$smCounter2_TP53_VAF[
  scatter_VAF$smCounter2_TP53_VAF != 0 & scatter_VAF$smCounter2_TP53_VAF != "not_detected" ])))
plot(hist(as.numeric(scatter_VAF$smCounter2_STAG2_VAF[
  scatter_VAF$smCounter2_STAG2_VAF != 0 & scatter_VAF$smCounter2_STAG2_VAF != "not_detected" ])))
dev.off()

# take the mean of any 2 value VAFs:
plot_VAFs <- scatter_VAF
plot_VAFs$GeneGlobe_TP53_VAF[grep("/", plot_VAFs$GeneGlobe_TP53_VAF)] <- 
  sapply(strsplit(
    plot_VAFs$GeneGlobe_TP53_VAF[grep("/", plot_VAFs$GeneGlobe_TP53_VAF)], "/" ),
    function(x) mean(as.numeric(x)) )
plot_VAFs$GeneGlobe_STAG2_VAF[grep("/", plot_VAFs$GeneGlobe_STAG2_VAF)] <- 
  sapply(strsplit(
    plot_VAFs$GeneGlobe_STAG2_VAF[grep("/", plot_VAFs$GeneGlobe_STAG2_VAF)], "/" ),
    function(x) mean(as.numeric(x)) )

# create smCounter2 TP53 vs STAG2 VAF plot:
VAF_df <- subset(scatter_VAF, select = c(
  Patient_id, Treatment.dilution, smCounter2_TP53_VAF, smCounter2_STAG2_VAF ))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment[VAF_df$id == "ES8"] <- "ES8"
VAF_df$treatment[VAF_df$id == "A673"] <- "A673"
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
VAF_df <- subset(scatter_VAF, select = c(
  Patient_id, Treatment.dilution, smCounter2_TP53_VAF, VAF_forward ))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment[VAF_df$id == "ES8"] <- "ES8"
VAF_df$treatment[VAF_df$id == "A673"] <- "A673"

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
VAF_df <- subset(scatter_VAF, select = c(
  Patient_id, Treatment.dilution, smCounter2_STAG2_VAF, VAF_forward ))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment[VAF_df$id == "ES8"] <- "ES8"
VAF_df$treatment[VAF_df$id == "A673"] <- "A673"

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

