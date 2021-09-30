
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

# fetch and annotate VAFs:
sm_vars <- list(
  TP53 = do.call("rbind", lapply(meta_list, fetch_sm_vafs, TP53)),
  STAG2 = do.call("rbind", lapply(meta_list, fetch_sm_vafs, STAG2)) )

# format:
sm_vars <- lapply(sm_vars, function(x))
for (i in seq_along(sm_vars)) {
  sm_vars[[i]]$new_SNV <- paste0(sm_vars[[i]]$seqnames, ":", sm_vars[[i]]$start, "_", 
    sm_vars[[i]]$ref, ">", sm_vars[[i]]$alt )
  var_df <- subset(sm_vars[[i]], select = c(new_SNV, effect_size, VAF) )
  colnames(var_df) <- c("smCounter2_SNV", "smCounter2_effect_size", "smCounter2_VAF")
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

# remove SNV info for not_detected samples:
final_vars$smCounter2_TP53_SNV[grep("not_detected", final_vars$smCounter2_TP53_SNV)] <- 
  "not_detected"
final_vars$smCounter2_STAG2_SNV[grep("not_detected", final_vars$smCounter2_STAG2_SNV)] <- 
  "not_detected"

# subset:
final_df <- subset(final_vars, 
  select = c(Patient_id, Sample_id, Library_id, Site, Treatment, 
    Sanger_TP53_SNV, TP53_SNV_type, smCounter2_TP53_SNV, ddPCR_TP53_VAF, 
    GeneGlobe_TP53_VAF, smCounter2_TP53_VAF, smCounter2_TP53_effect_size,
    Sanger_STAG2_SNV, smCounter2_STAG2_SNV, ddPCR_STAG2_VAF, 
    GeneGlobe_STAG2_VAF, smCounter2_STAG2_VAF, smCounter2_STAG2_effect_size,
    Pathology_EWSR1_FLI1 ) )

# sort by original metadata order:
final_df <- final_df[match(meta$Sample_id, final_df$Sample_id),]

write.table(
  final_df,
  file.path(table_dir, "sample_summary_with_smcounter_SNV.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE )



####################################################################################
### 1. Plot GeneGlobe vs smcounter2 SNP VAFs ###
####################################################################################

# plot GG_TP53 vs SM_TP53 VAFs, leaving out samples with unknown results:
GG_vs_SM_TP53 <- subset(
  patient_VAFs, select = c(GG_TP53_VAF, SM_TP53_VAF, Treatment, Patient)
)
GG_vs_SM_TP53$Type <- "treated ctDNA"
GG_vs_SM_TP53$Type[GG_vs_SM_TP53$Treatment == "tumour"] <- "gDNA"
GG_vs_SM_TP53$Type[GG_vs_SM_TP53$Treatment == "naive"] <- "untreated ctDNA"
GG_vs_SM_TP53 <- GG_vs_SM_TP53[
  GG_vs_SM_TP53$GG_TP53_VAF != "unknown" & 
  GG_vs_SM_TP53$SM_TP53_VAF != "unknown",
]

# convert numbers to numeric:
GG_vs_SM_TP53$GG_TP53_VAF <- as.numeric(GG_vs_SM_TP53$GG_TP53_VAF)
GG_vs_SM_TP53$SM_TP53_VAF <- as.numeric(GG_vs_SM_TP53$SM_TP53_VAF)

p <- ggplot(
  GG_vs_SM_TP53, 
  aes(x = GG_TP53_VAF, y = SM_TP53_VAF, color = Treatment)
)
p <- p + geom_point(size = 3)
p <- p + theme_cowplot(12)
p <- p + ylim(c(0, 100))
p <- p + xlim(c(0, 100))
p <- p + ylab("smCounter2 TP53 VAF")
p <- p + xlab("GeneGlobe TP53 VAF")
p <- p + geom_text_repel(data=GG_vs_SM_TP53, aes(label=Patient))

png(
  paste0(plot_dir, "GG_TP53_vs_SM_TP53_VAFs.png"),
  height = 3.5,
  width = 5,
  res = 300,
  units = "in")
  print(p)
dev.off()

# plot GG_STAG2 vs SM_STAG2 VAFs, leaving out samples with unknown results:
GG_vs_SM_STAG2 <- subset(
  patient_VAFs, select = c(GG_STAG2_VAF, SM_STAG2_VAF, Treatment, Patient)
)
GG_vs_SM_STAG2$Type <- "treated ctDNA"
GG_vs_SM_STAG2$Type[GG_vs_SM_STAG2$Treatment == "tumour"] <- "gDNA"
GG_vs_SM_STAG2$Type[GG_vs_SM_STAG2$Treatment == "naive"] <- "untreated ctDNA"
GG_vs_SM_STAG2 <- GG_vs_SM_STAG2[
  GG_vs_SM_STAG2$GG_STAG2_VAF != "unknown" & 
    GG_vs_SM_STAG2$SM_STAG2_VAF != "unknown",
]

# convert numbers to numeric:
GG_vs_SM_STAG2$GG_STAG2_VAF <- as.numeric(GG_vs_SM_STAG2$GG_STAG2_VAF)
GG_vs_SM_STAG2$SM_STAG2_VAF <- as.numeric(GG_vs_SM_STAG2$SM_STAG2_VAF)

p <- ggplot(
  GG_vs_SM_STAG2, 
  aes(x = GG_STAG2_VAF, y = SM_STAG2_VAF, color = Treatment)
)
p <- p + geom_point(size = 3)
p <- p + theme_cowplot(12)
p <- p + ylim(c(0, 100))
p <- p + xlim(c(0, 100))
p <- p + ylab("smCounter2 STAG2 VAF")
p <- p + xlab("GeneGlobe STAG2 VAF")
p <- p + geom_text_repel(data=GG_vs_SM_STAG2, aes(label=Patient))

png(
  paste0(plot_dir, "GG_STAG2_vs_SM_STAG2_VAFs.png"),
  height = 3.5,
  width = 5,
  res = 300,
  units = "in")
print(p)
dev.off()


####################################################################################
### 2. Plot TP53 vs STAG2 VAFs ###
####################################################################################

# plot TP53 vs STAG2 GG VAFs, leaving out samples with either = 0 or "unknown"
GG_TP53_vs_STAG2 <- subset(
  patient_VAFs, select = c(GG_TP53_VAF, GG_STAG2_VAF, Treatment, Patient)
)
GG_TP53_vs_STAG2 <- GG_TP53_vs_STAG2[
  GG_TP53_vs_STAG2$GG_TP53_VAF != 0 & 
  GG_TP53_vs_STAG2$GG_STAG2_VAF != 0 &
  GG_TP53_vs_STAG2$GG_TP53_VAF != "unknown" & 
  GG_TP53_vs_STAG2$GG_STAG2_VAF != "unknown",
]

# convert to numeric values:
GG_TP53_vs_STAG2$GG_TP53_VAF <- as.numeric(
  GG_TP53_vs_STAG2$GG_TP53_VAF
)
GG_TP53_vs_STAG2$GG_STAG2_VAF <- as.numeric(
  GG_TP53_vs_STAG2$GG_STAG2_VAF
)

p <- ggplot(
  GG_TP53_vs_STAG2, 
  aes(x = GG_TP53_VAF, y = GG_STAG2_VAF, color = Treatment)
)
p <- p + geom_point(size = 3)
p <- p + theme_cowplot(12)
p <- p + ylim(c(0, 100))
p <- p + xlim(c(0, 100))
p <- p + ylab("GeneGlobe STAG2 VAF")
p <- p + xlab("GeneGlobe TP53 VAF")
p <- p + geom_text_repel(data=GG_TP53_vs_STAG2, aes(label=Patient))

png(
  paste0(plot_dir, "GG_TP53_vs_STAG2_VAFs.png"),
  height = 3.5,
  width = 5,
  res = 300,
  units = "in")
  print(p)
dev.off()

# plot TP53 vs STAG2 SM VAFs, leaving out samples with either = 0
SM_TP53_vs_STAG2 <- subset(
  patient_VAFs, select = c(SM_TP53_VAF, SM_STAG2_VAF, Treatment, Patient)
)

SM_TP53_vs_STAG2 <- SM_TP53_vs_STAG2[
  SM_TP53_vs_STAG2$SM_TP53_VAF != 0 & 
    SM_TP53_vs_STAG2$SM_STAG2_VAF != 0 &
    SM_TP53_vs_STAG2$SM_TP53_VAF != "unknown" & 
    SM_TP53_vs_STAG2$SM_STAG2_VAF != "unknown",
]

# calculate correlation between TP53 and STAG2 VAFs:
plot(hist(SM_TP53_vs_STAG2$SM_TP53_VAF))
plot(hist(SM_TP53_vs_STAG2$SM_STAG2_VAF))
dev.off()

corr <- cor.test(
  x = SM_TP53_vs_STAG2$SM_TP53_VAF, 
  y = SM_TP53_vs_STAG2$SM_STAG2_VAF,
  method = "pearson"
)

# Fit regression line
require(stats)
reg <- lm(SM_TP53_VAF ~ SM_STAG2_VAF, data = SM_TP53_vs_STAG2)
coeff=coefficients(reg)

p <- ggplot(
  SM_TP53_vs_STAG2, 
  aes(x = SM_TP53_VAF, y = SM_STAG2_VAF, color = Treatment)
)
p <- p + geom_point(size = 3)
p <- p + theme_cowplot(12)
p <- p + ylim(c(0, 100))
p <- p + xlim(c(0, 100))
p <- p + ylab("smCounter2 STAG2 VAF")
p <- p + xlab("TP53 VAF")
p <- p + geom_text_repel(data=SM_TP53_vs_STAG2, aes(label=Patient))
p <- p + annotate(
  "text", x = 55, y = 90, 
  label = paste0(
    "R2=", round(corr$estimate, 2), ", p=", round(corr$p.value, 3)
  ), color='red', size = 4
)
p <- p + geom_abline(intercept = coeff[1], slope = coeff[2], color = "red")

png(paste0(plot_dir, "SM_TP53_vs_STAG2_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(p)
dev.off()


####################################################################################
### 2. Plot fusion vs GeneGlobe SNP VAFs ###
####################################################################################

# plot fusion vs TP53 VAFs, leaving out samples with either = 0 or unknown:
fusion_vs_SM_TP53 <- subset(
  patient_VAFs, select = c(SM_TP53_VAF, fusion_VAF, Treatment, Patient)
)

fusion_vs_SM_TP53 <- fusion_vs_SM_TP53[
  fusion_vs_SM_TP53$SM_TP53_VAF != 0 & 
    fusion_vs_SM_TP53$fusion_VAF != 0 &
    fusion_vs_SM_TP53$SM_TP53_VAF != "unknown" & 
    fusion_vs_SM_TP53$fusion_VAF != "unknown",
]

# calculate correlation between fusion and TP53 VAFs:
plot(hist(fusion_vs_SM_TP53$fusion_VAF))
plot(hist(fusion_vs_SM_TP53$SM_TP53_VAF))
dev.off()

corr <- cor.test(
  x = fusion_vs_SM_TP53$SM_TP53_VAF, 
  y = fusion_vs_SM_TP53$fusion_VAF,
  method = "pearson"
)

# Fit regression line
require(stats)
reg <- lm(SM_TP53_VAF ~ fusion_VAF, data = fusion_vs_SM_TP53)
coeff=coefficients(reg)

p <- ggplot(
  fusion_vs_SM_TP53, 
  aes(x = fusion_VAF, y = SM_TP53_VAF, color = Treatment)
)
p <- p + geom_point(size = 3)
p <- p + theme_cowplot(12)
p <- p + ylim(c(0, 100))
p <- p + xlim(c(0, 100))
p <- p + ylab("smCounter2 TP53 VAF")
p <- p + xlab("fusion VAF")
p <- p + geom_text_repel(data=fusion_vs_SM_TP53, aes(label=Patient))
p <- p + annotate(
  "text", x = 80, y = 85, 
  label = paste0(
    "R2=", round(corr$estimate, 2), ", p=", round(corr$p.value, 3)
  ), color='red', size = 4
)
p <- p + geom_abline(intercept = coeff[1], slope = coeff[2], color = "red")

png(paste0(plot_dir, "fusion_vs_SM_TP53_VAFs.png"),
  height = 3.5,
  width = 5,
  res = 300,
  units = "in")
  print(p)
dev.off()

# plot fusion vs STAG2 VAFs, leaving out samples with either = 0 or unknown:
fusion_vs_SM_STAG2 <- subset(
  patient_VAFs, select = c(SM_STAG2_VAF, fusion_VAF, Treatment, Patient)
)

fusion_vs_SM_STAG2 <- fusion_vs_SM_STAG2[
  fusion_vs_SM_STAG2$SM_STAG2_VAF != 0 & 
    fusion_vs_SM_STAG2$fusion_VAF != 0 &
    fusion_vs_SM_STAG2$SM_STAG2_VAF != "unknown" & 
    fusion_vs_SM_STAG2$fusion_VAF != "unknown",
]
# calculate correlation between fusion and STAG2 VAFs:
plot(hist(fusion_vs_SM_STAG2$fusion_VAF))
plot(hist(fusion_vs_SM_STAG2$SM_STAG2_VAF))
dev.off()

corr <- cor.test(
  x = fusion_vs_SM_STAG2$SM_STAG2_VAF, 
  y = fusion_vs_SM_STAG2$fusion_VAF,
  method = "pearson"
)

# Fit regression line
require(stats)
reg <- lm(SM_STAG2_VAF ~ fusion_VAF, data = fusion_vs_SM_STAG2)
coeff=coefficients(reg)

p <- ggplot(
  fusion_vs_SM_STAG2, 
  aes(x = fusion_VAF, y = SM_STAG2_VAF, color = Treatment)
)
p <- p + geom_point(size = 3)
p <- p + theme_cowplot(12)
p <- p + ylim(c(0, 100))
p <- p + xlim(c(0, 100))
p <- p + ylab("smCounter2 STAG2 VAF")
p <- p + xlab("fusion VAF")
p <- p + geom_text_repel(data=fusion_vs_SM_STAG2, aes(label=Patient))
p <- p + annotate(
  "text", x = 80, y = 85, 
  label = paste0(
    "R2=", round(corr$estimate, 2), ", p=", round(corr$p.value, 3)
  ), color='red', size = 4
)
p <- p + geom_abline(intercept = coeff[1], slope = coeff[2], color = "red")

png(paste0(plot_dir, "fusion_vs_SM_STAG2_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(p)
dev.off()


####################################################################################
### 3. Plot fusion vs smcounter2 SNP VAFs ###
####################################################################################


