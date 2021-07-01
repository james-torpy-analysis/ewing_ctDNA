
home_dir <- "/Users/torpor/clusterHome/"
#home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
result_dir <- paste0(project_dir, "results/")
in_path <- paste0(result_dir, "VAF_calculation/")


####################################################################################
### 0. Load packages and functions ###
####################################################################################

library(ggplot2)
library(cowplot)


####################################################################################
### 1. Load VAFs ###
####################################################################################

samplenames <- as.list(list.files(in_path, pattern = "409"))

VAFs <- sapply(samplenames, function(x) {
  return(
    read.table(
      paste0(in_path, x, "/VAF.txt"),
      sep = "\t",
      header = T,
      stringsAsFactors = F
    )
  )
})

######
VAF_df <- data.frame(
  id = c(
    "409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001", 
    "409_032_DCB94_GGACTCCT-CTCTCTAT_L001",
    "409_033_DCB94_TAGGCATG-CTCTCTAT_L001"
  ),
  dilution = factor(c(100, 50, 10), levels = c(100, 50, 10)),
  VAF = c(40.4, 30.4, 20.4),
  fusion_coord = c(29684299, 29684300, 29684301)
)
######


####################################################################################
### 2. Plot VAFs ###
####################################################################################

p <- ggplot(VAF_df, aes(x = dilution, y = VAF))
p <- p + geom_bar(stat="identity")
p <- p + scale_x_discrete()
p <- p + theme_cowplot(12)
p <- p + xlab("dilution (%)")
p <- p + ylab("VAF (%)")
p <- p + geom_text(
  aes(label=VAF), 
  vjust=-0.25,
  fontface = "bold"
)
p

png(paste0(in_path, "serial_dilution_1_VAFs.png"))
  p
dev.off()

