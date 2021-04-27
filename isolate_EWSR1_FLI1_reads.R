#!/bin/bash

sample_name <- "409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
in_path <- paste0(project_dir, "results/picard/bams/", sample_name, "/")

func_dir <- paste0(project_dir, "scripts/functions/")

# define GOI co-ordinates:
GOI <- list(
  EWSR1 = GRanges(
    seqnames = "chr22",
    ranges = IRanges(start = 29664257, end = 29696511),
    strand = "*"
  ),
  FLI1 = GRanges(
    seqnames = "chr11",
    ranges = IRanges(start = 128554430, end = 128685162),
    strand = "*"
  )
)


####################################################################################
### 0. Load packages ###
####################################################################################

library(rtracklayer)
library(Rsamtools)


####################################################################################
### 1. Isolate reads across chr11 and chr22 ###
####################################################################################

system(
  paste0(
    "samtools view ", in_path, sample_name, 
    ".uncollapsed.bam | grep chr11 | grep chr22 > ", 
    in_path, sample_name, ".chr11_22.uncollapsed.sam"
  )
)


####################################################################################
### 1. Load bam ###
####################################################################################

#param = ScanBamParam(
#  what = c(
#    "qname", "rname", "pos"
#  )
#)  
#bam <- scanBam(
#  paste0(in_path, sample_name, ".chr11_22.uncollapsed.bam"),
#  param = param
#)
#
## remove NAs:
#not_na <- which(!is.na(bam[[1]]$pos))
#filt_bam <- lapply(bam[[1]], function(x) x[not_na])
#
#bam_gr <- GRanges(
#  seqnames = filt_bam$qname,
#  ranges = IRanges(start = filt_bam$pos, end = filt_bam$pos),
#  strand = "*",
#  rname = filt_bam$rname
#)

bam <- read.table(
  paste0(in_path, sample_name, ".chr11_22.uncollapsed.sam"),
  sep = "\t",
  fill = TRUE,
  header = F,
  stringsAsFactors = F
)

bam_gr <- GRanges(
  seqnames = bam$V3,
  ranges = IRanges(start = bam$V4, end = bam$V4),
  strand = "*",
  rname = bam$V1,
  flag = bam$V2
)

bam_mate_gr <- GRanges(
  seqnames = bam$V7,
  ranges = IRanges(start = bam$V8, end = bam$V8),
  strand = "*",
  rname = bam$V1,
  flag = bam$V2
)


####################################################################################
### 2. Find overlaps of EWSR1 and FLI1 ###
####################################################################################

EWSR1_olaps <- findOverlaps(bam_gr, GOI$EWSR1)
EWSR1_gr <- bam_gr[queryHits(EWSR1_olaps)]

EWSR1_mate_olaps <- findOverlaps(bam_mate_gr, GOI$EWSR1)
EWSR1_mate_gr <- bam_mate_gr[queryHits(EWSR1_mate_olaps)]

FLI1_olaps <- findOverlaps(bam_gr, GOI$FLI1)
FLI1_gr <- bam_gr[queryHits(FLI1_olaps)]

FLI1_mate_olaps <- findOverlaps(bam_mate_gr, GOI$FLI1)
FLI1_mate_gr <- bam_mate_gr[queryHits(FLI1_mate_olaps)]

fusions <- c(
  EWSR1_gr$rname[EWSR1_gr$rname %in% FLI1_mate_gr$rname],
  EWSR1_mate_gr$rname[EWSR1_mate_gr$rname %in% FLI1_gr$rname]
)


####################################################################################
### 3. Isolate fusion reads from bam ###
####################################################################################

# insert header into new bam:
system(
  paste0(
    "samtools view -H ", in_path, sample_name, ".uncollapsed.bam > ",
    in_path, sample_name, ".EWSR1_FLI1_fusion.uncollapsed.sam"
  )
)

# grep each fusion read name to put into output file:
for (rname in fusions) {
  system(
    paste0(
      "grep ", rname, " ", in_path, sample_name, ".chr11_22.uncollapsed.sam", 
      " >> ", in_path, sample_name, ".EWSR1_FLI1_fusion.uncollapsed.sam"
    )
  )
}

# convert to bam:
system(
  paste0(
    "samtools view -bh ", in_path, sample_name, ".EWSR1_FLI1_fusion.uncollapsed.sam > ",
    in_path, sample_name, ".EWSR1_FLI1_fusion.uncollapsed.bam"
  )
)

# sort and index:
system(
  paste0(
    "samtools sort -o ", in_path, sample_name, ".EWSR1_FLI1_fusion.uncollapsed.sorted.bam ", 
    in_path, sample_name, ".EWSR1_FLI1_fusion.uncollapsed.bam"
  )
)

system(
  paste0(
    "samtools index ", in_path, sample_name, ".EWSR1_FLI1_fusion.uncollapsed.sorted.bam"
  )
)

# clean up sams:
system(
  paste0(
    "rm ", in_path, sample_name, ".chr11_22.uncollapsed.sam ", 
    in_path, sample_name, ".EWSR1_FLI1_fusion.uncollapsed.sam"
  )
)



