
summary_file_order <- c(
  "001", "002", "050", "008", "010", 
  "003", "052", "021", "024", "026", 
  "028", "055", "004", "051", "005", 
  "056", "009", "011", "006", "007", 
  "012", "013", "017", "049", "014", 
  "060", "015", "016", "054", "018", 
  "023", "027", "019", "059", "048", 
  "020", "057", "022", "030", "031", 
  "062", "025", "058", "061", "063", 
  "002", "040", "041", "065", "066", 
  "067", "042", "043", "044", "045", 
  "046", "047", "001", "032", "033", 
  "068", "069", "070", "034", "035", 
  "036", "037", "038", "039"
)

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")

detection_dir <- paste0(project_dir, "results/detection_heatmaps/")
VAF_dir <- paste0(detection_dir, "/Rdata/")
variant_path <- paste0(project_dir, "results/smcounter2/")

plot_dir <- paste0(detection_dir, "plots/")
table_dir <- paste0(detection_dir, "tables/")
Robject_dir <- paste0(detection_dir, "Rdata/")

system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", Robject_dir))

library(ggplot2)
library(GenomicAlignments)
library(tibble)
library(cowplot)
library(ggrepel)

fetch_sm_vafs <- dget(paste0(func_dir, "fetch_sm_vafs.R"))


####################################################################################
### 0. Load mutation data ###
####################################################################################

# load fusion VAFs:
fusion_VAFs <- readRDS(paste0(VAF_dir, "all_VAFs.Rdata"))
rownames(fusion_VAFs) <- fusion_VAFs$Sample
fusion_VAFs$VAF <- round(fusion_VAFs$VAF*100, 1)
colnames(fusion_VAFs) <- "fusion_VAF"

# load and add SNP VAFs:
patient_meta <- read.table(
  paste0(ref_dir, "/ES_samples_by_patient.tsv"),
  header = T,
  fill = T
)
# label SNP VAFs as GeneGlobe:
colnames(patient_meta) <- gsub("TP53", "GG_TP53", colnames(patient_meta))
colnames(patient_meta) <- gsub("STAG2", "GG_STAG2", colnames(patient_meta))
patient_VAFs <- subset(
  patient_meta, select = c(
    "Patient", "Sample", "Treatment", "GG_TP53_VAF", "GG_STAG2_VAF"
  )
)
patient_VAFs$GG_TP53_VAF <- gsub("/.*$", "", patient_VAFs$GG_TP53_VAF)
rownames(patient_VAFs) <- gsub("_.*$", "", 
  gsub("409_", "", patient_VAFs$Sample))
patient_VAFs <- merge(patient_VAFs, fusion_VAFs, by=0, all = FALSE)
patient_VAFs <- patient_VAFs[,!is.na(colnames(patient_VAFs))]

## load smcounter2 SNP VAFs:
# define samplenames:
samplenames <- as.list(list.files(variant_path, pattern = "409"))

# define TP53 and STAG2 ROI:
snp_roi <- read.table(
  paste0(ref_dir, "CDHS-34925Z-409.roi.bed"),
  sep = "\t"
)
roi_gr <- GRanges(
  seqnames = snp_roi$V1,
  ranges = IRanges(start = snp_roi$V2, end = snp_roi$V3),
  gene = snp_roi$V4
)

# split roi into genes and reduce to one set coords:
roi_spl <- split(roi_gr, roi_gr$gene)
roi <- lapply(roi_spl, function(x) {
  GRanges(
    seqnames = seqnames(x)[1],
    ranges = IRanges(start = start(x)[1], end = end(x)[length(x)]),
    strand = "*",
    gene = x$gene[1]
  )
})

# add Treatment information to samplenames:
sample_dfs <- split(subset(patient_VAFs, select = c(Sample, Treatment)), 
  patient_VAFs$Sample )

# fetch and annotate VAFs:
sm_vars <- lapply(sample_dfs, fetch_sm_vafs, roi, alt_var)
names(sm_vars) <- patient_VAFs$Sample

# fetch VAF vectors:
TP53_VAF <- lapply(sm_vars, function(x) {
  data.frame(SM_TP53_VAF = x$TP53$VAF, effect_size = x$TP53$effect_size) })
TP53_VAF <- do.call("rbind", TP53_VAF)
rownames(TP53_VAF) <- gsub(
  "_.*$", "", gsub("409_", "", rownames(TP53_VAF))
)

STAG2_VAF <- unlist(lapply(sm_vars, function(x) x$STAG2$VAF))
STAG2_VAF <- data.frame(
  row.names = gsub(
    "_.*$", "", gsub("409_", "", names(STAG2_VAF))
  ),
  SM_STAG2_VAF = STAG2_VAF
)

save.image(paste0(Robject_dir, "adding_VAFs.Rdata"))

# merge with patient VAFs:
patient_VAFs <- patient_VAFs %>%
  column_to_rownames("Row.names")
patient_VAFs <- merge(patient_VAFs, TP53_VAF, by=0, all = TRUE)

patient_VAFs <- patient_VAFs %>%
  column_to_rownames("Row.names")
patient_VAFs <- merge(patient_VAFs, STAG2_VAF, by=0, all = TRUE)

patient_VAFs <- patient_VAFs[!(is.na(patient_VAFs$Patient)),]

# make NA values = 0:
patient_VAFs$SM_TP53_VAF[is.na(patient_VAFs$SM_TP53_VAF)] <- 0
patient_VAFs$SM_STAG2_VAF[is.na(patient_VAFs$SM_STAG2_VAF)] <- 0

# write table:
colnames(patient_VAFs)[1] <- "ID"
write.table(
  patient_VAFs, 
  paste0(table_dir, "patient_VAFs.tsv"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T)

# order for summary file:
summary_file_order <- summary_file_order[
  summary_file_order %in% patient_VAFs$ID]
patient_VAFs <- patient_VAFs[match(summary_file_order, patient_VAFs$ID),]

# write table:
write.table(
  patient_VAFs, 
  paste0(table_dir, "patient_VAFs_summary_order.tsv"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T)


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


