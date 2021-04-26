#!/bin/bash

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
in_path <- paste0(project_dir, "results/svaba/bwa/")

func_dir <- paste0(project_dir, "scripts/functions/")

# define GOI co-ordinates:
GOI <- list(
  FLI1 = GRanges(
    seqnames = "chr11",
    ranges = IRanges(start = 128554430, end = 128685162),
    strand = "*"
  ),
  ERG = GRanges(
    seqnames = "chr7",
    ranges = IRanges(start = 13928854, end = 14033050),
    strand = "*"
  ),
  ETV1 = GRanges(
    seqnames = "chr21",
    ranges = IRanges(start = 39737183, end = 40035704),
    strand = "*"
  ),
  ETV4 = GRanges(
    seqnames = "chr17",
    ranges = IRanges(start = 41603214, end = 41625708),
    strand = "*"
  ),
  FEV = GRanges(
    seqnames = "chr2",
    ranges = IRanges(start = 219843809, end = 219851906),
    strand = "*"
  )
)


####################################################################################
### 0. Load packages ###
####################################################################################

library(rtracklayer)
library(VariantAnnotation)

find_fusions <- dget(paste0(func_dir, "find_fusions.R"))


####################################################################################
### 1. Load VCFs and define EWSR1 co-ordinates ###
####################################################################################

# fetch sample names:
samplenames <- grep(
  "_001_|_002_",
  list.files(in_path, pattern = "409"),
  invert = T,
  value = T
)

# load VCFs and convert to granges:
breakpoint_gr <- lapply(samplenames, function(x) {

  print(paste0("Loading ", x, "..."))

  vcf_df <- tryCatch(
  	read.table(
      paste0(in_path, x, "/", x, ".svaba.unfiltered.sv.formatted.vcf"),
  	  sep = "\t"
    ),
    warning = function(w) {NA},
    error = function(e) {NA}
  )
  
  if (!is.na(vcf_df)) {

  	# split additional info field:
  	additional_info <- gsub("^.*=", "", do.call("rbind", strsplit(vcf_df$V8, ";")))

    return(
      GRanges(
        seqnames = vcf_df$V1,
        ranges = IRanges(start = vcf_df$V2, vcf_df$V2),
        strand = "*",
        id = vcf_df$V3,
        join_base = vcf_df$V4,
        join = vcf_df$V5,
        quality = vcf_df$V6,
        filter = vcf_df$V7,
        DISC_MAPQ = additional_info[,1],
        EVDNC = additional_info[,2],
        type = additional_info[,3],
        MAPQ = additional_info[,4],
        MATEID = additional_info[,5],
        MATENM = additional_info[,6],
        NM = additional_info[,7],
        NUMPARTS = additional_info[,8]
      )
    )

  } else {
    return(NA)
  }

})
names(breakpoint_gr) <- samplenames


####################################################################################
### 2. Find overlaps with breakpoints and EWSR1 ###
####################################################################################

# define EWSR1 co-ords:
EWSR1 <- GRanges(
  seqnames = "chr22",
  ranges = IRanges(start = 29663998, end = 29696515),
  strand = "*"
)

# find breakpoint overlaps with EWSR1:
EWSR1_fusions <- lapply(breakpoint_gr, find_fusions, EWSR1)

# remove NAs:
EWSR1_fusions <- EWSR1_fusions[!is.na(EWSR1_fusions)]

# count ranges:
print("EWSR1 fusions per sample:")
print(sapply(EWSR1_fusions, length))

# determine whether joining ranges overlap FLI1, ETV1 or ERG:
for (i in 1:length(EWSR1_fusions)) {

  print(i)

  if (length(EWSR1_fusions[[i]]) > 0) {

    seqnams <- gsub(
      ":.*$", "", 
      gsub("^.*chr", "chr", EWSR1_fusions[[i]]$join)
    )
    coord <- as.numeric(
      gsub(
        "[^0-9.-]", "", 
        gsub("^.*chr.*:", "", EWSR1_fusions[[i]]$join)
      )
    )
  
    if (!exists("join_ranges")) {
      join_ranges <- list(
        GRanges(
          seqnames = seqnams,
          ranges = IRanges(start = coord, end = coord),
          strand = "*",
          join_chr = "chr22",
          join = start(EWSR1_fusions[[i]])
        )
      )
    } else {
      join_ranges[[i]] <- GRanges(
        seqnames = seqnams,
        ranges = IRanges(start = coord, end = coord),
        strand = "*",
        join_chr = "chr22",
        join_coord = start(EWSR1_fusions[[i]])
      )
    }

  } else {

  	if (!exists("join_ranges")) {
      join_ranges <- list(NA)
    } else {
      join_ranges[[i]] <- NA
    }

  }  

}
names(join_ranges) <- names(EWSR1_fusions)

# identify EWSR1 fusions with GOI: 
GOI_fusions <- lapply(GOI, function(x) {

  for (i in 1:length(join_ranges)) {
  
    if (i==1) {
      fusions <- list(find_fusions(join_ranges[[i]], x))
    } else {
      fusions[[i]] <- find_fusions(join_ranges[[i]], x)
    }
  
  }
  names(fusions) <- names(join_ranges)

  return(fusions)

})

# count GOI fusions:
fusion_nos <- lapply(GOI_fusions, function(x) {

  lengths <- sapply(x, function(y) {
    if (!is.na(y)) {
      return(length(y))
    } else {
      return(0)
    }
  })
  names(lengths) <- gsub(
    "_.*$", "", 
    gsub("409_", "", names(x))
  )

  return(lengths)

})
