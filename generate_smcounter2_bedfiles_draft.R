home_dir <- "/share/ScratchGeneral/jamtor"
project_dir <- file.path(home_dir, "projects/ewing_ctDNA")
in_dir <- file.path(project_dir, "refs")

library(GenomicRanges)

# load primers:
primers <- read.table(
  file.path(in_dir, "CDHS-34925Z-409.primers.txt"),
  sep = "\t",
  header = FALSE )

# annotate orientation:
primers$V3[primers$V3 == 0] <- "-"
primers$V3[primers$V3 == 1] <- "+"

# convert to gr:
gr <- GRanges(
  seqnames <- primers$V1,
  ranges <- IRanges(start = primers$V2, width = length(primers$V4)),
  strand <- primers$V3 )
