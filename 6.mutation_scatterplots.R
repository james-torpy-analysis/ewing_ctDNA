
#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
ref_dir <- paste0(project_dir, "refs/")

detection_dir <- paste0(project_dir, "results/detection_heatmaps/")
VAF_dir <- paste0(detection_dir, "/Rdata/")
out_dir <- paste0(detection_dir, "plots/")

system(paste0("mkdir -p ", out_dir))

library(ggplot2)

####################################################################################
### 0. Load mutation data ###
####################################################################################

# load fusion VAFs:
fusion_VAFs <- readRDS(paste0(VAF_dir, "all_VAFs.Rdata"))
fusion_VAFs$VAF <- round(fusion_VAFs$VAF*100, 1)
colnames(fusion_VAFs) <- "fusion_VAF"

# load and add SNP VAFs:
patient_meta <- read.table(
  paste0(ref_dir, "/ES_samples_by_patient.txt"),
  header = T,
  fill = T
)
patient_VAFs <- subset(
  patient_meta, select = c(
    "Patient", "Sample", "Treatment", "GG_TP53_VAF", "GG_STAG2_VAF"
  )
)
patient_VAFs$GG_TP53_VAF <- gsub("/.*$", "", patient_VAFs$GG_TP53_VAF)
rownames(patient_VAFs) <- gsub("_.*$", "", 
  gsub("409_", "", patient_VAFs$Sample))
patient_VAFs <- merge(patient_VAFs, fusion_VAFs, by=0, all = FALSE)


####################################################################################
### 1. Plot VAFs vs GeneGlobe variant ###
####################################################################################


# convert numbers to numeric:
patient_VAFs$GG_TP53_VAF <- as.numeric(patient_VAFs$GG_TP53_VAF)
patient_VAFs$GG_STAG2_VAF <- as.numeric(patient_VAFs$GG_STAG2_VAF)
patient_VAFs$fusion_VAF <- as.numeric(patient_VAFs$fusion_VAF)

# plot TP53 vs STAG2 VAFs, leaving out samples with either = 0
patient_TP53_vs_STAG2 <- subset(
  patient_VAFs, select = c(GG_TP53_VAF, GG_STAG2_VAF, Treatment)
)
patient_TP53_vs_STAG2 <- patient_TP53_vs_STAG2[
  patient_TP53_vs_STAG2$GG_TP53_VAF != 0 & patient_TP53_vs_STAG2$GG_STAG2_VAF != 0,
]


p <- ggplot(
  patient_TP53_vs_STAG2, aes(x = GG_TP53_VAF, y = GG_STAG2_VAF, color = Treatment)
)
p <- p + geom_point()

png(paste0(out_dir, "TP53_vs_GG_STAG2_VAFs.png"))
  print(p)
dev.off()

# plot fusion vs TP53 VAFs, leaving out samples with either = 0
patient_fusion_vs_TP53 <- subset(
  patient_VAFs, select = c(GG_TP53_VAF, fusion_VAF, Treatment)
)
patient_fusion_vs_TP53 <- patient_fusion_vs_TP53[
  patient_fusion_vs_TP53$GG_TP53_VAF != 0 & patient_fusion_vs_TP53$fusion_VAF != 0,
]

p <- ggplot(
  patient_fusion_vs_TP53, aes(x = fusion_VAF, y = GG_TP53_VAF, color = Treatment)
)
p <- p + geom_point()
png(paste0(out_dir, "fusion_vs_GG_TP53_VAFs.png"))
  print(p)
dev.off()

# plot fusion vs STAG2 VAFs, leaving out samples with either = 0
patient_fusion_vs_STAG2 <- subset(
  patient_VAFs, select = c(GG_STAG2_VAF, fusion_VAF, Treatment)
)
patient_fusion_vs_STAG2 <- patient_fusion_vs_STAG2[
  patient_fusion_vs_STAG2$GG_STAG2_VAF != 0 & patient_fusion_vs_STAG2$fusion_VAF != 0,
]

p <- ggplot(
  patient_fusion_vs_STAG2, aes(x = fusion_VAF, y = GG_STAG2_VAF, color = Treatment)
)
p <- p + geom_point()
png(paste0(out_dir, "fusion_vs_GG_STAG2_VAFs.png"))
  print(p)
dev.off()