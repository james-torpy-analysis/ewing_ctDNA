#!/bin/bash

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
in_path <- paste0(project_dir, "results/svaba/bwa/")


####################################################################################
### 0. Load packages ###
####################################################################################

library(rtracklayer)
library(VariantAnnotation)


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
#breakpoint_gr <- lapply(samplenames, function(x) {
#
#  vcf_df <- read.table(
#    paste0(in_path, x, "/", x, ".svaba.unfiltered.sv.formatted.vcf"),
#  	sep = "\t"
#  )
#  
#  return(
#  	vcf_gr <- GRanges(
#  	  seqnames = vcf_df$V1,
#  	  ranges = IRanges(start = vcf_df$V2, vcf_df$V2),
#  	  strand = "*",
#  	  id = vcf_df$V3,
#  	  join_base = vcf_df$V4,
#  	  join = vcf_df$V5,
#  	  quality = vcf_df$V6,
#  	  filter = vcf_df$V7,
#  	  additional_info = vcf_df$V8
#  	)
#  )
#
#})

for (s in 1:length(samplenames)) {

  print(samplenames[s])

  vcf_df <- tryCatch(
  	read.table(
      paste0(in_path, samplenames[s], "/", samplenames[s], ".svaba.unfiltered.sv.formatted.vcf"),
  	  sep = "\t"
    ),
    warning = function(w) {NA},
    error = function(e) {NA}
  )
  
  if (!is.na(vcf_df)) {

  	# split additional info field:
  	additional_info <- gsub("^.*=", "", do.call("rbind", strsplit(vcf_df$V8, ";")))

  	if (s==1) {
      vcf_gr <- list(
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
      vcf_gr[[s]] <- GRanges(
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
    }

  } else {

  	if (s==1) {
      vcf_gr <- list(
        NA
      )
    } else {
      vcf_gr[[s]] <- NA
    }

  }
  names(vcf_gr)[s] <- samplenames[s]

}

# define EWSR1 co-ords:
EWSR1 <- GRanges(
  seqnames = "chr22",
  ranges = IRanges(start = 29663998, end = 29696515),
  strand = "*"
)


####################################################################################
### 2. Find overlaps with breakpoints and EWSR1 ###
####################################################################################

######

find_fusions <- function(query_coord, subject_coord) {

  if (!is.na(query_coord)) {

    olaps <- findOverlaps(query_coord, subject_coord)

    if (length(olaps) > 0) {
      fusion <- query_coord[queryHits(olaps)]
    } else {
      fusion <- NA
    }

  } else {
    fusion <- NA
  }

  if (!is.na(fusion)) {
    # remove duplicates:
    fusion <- fusion[!duplicated(start(fusion))]
  
    # remove chr22 to chr22 fusions:
    fusion <- fusion[!( length(grep("chr22", seqnames(fusion))) > 0 & length(grep("chr22", fusion$join)) > 0 )]
  }

  return(fusion)
  
}

######

for (i in 1:length(vcf_gr)) {

  print(i)

  if (!is.na(vcf_gr[[i]])) {

  	olaps <- findOverlaps(vcf_gr[[i]], EWSR1)

    if (length(olaps) > 0) {
    	if (i==1) {
        EWSR1_fusions <- list(vcf_gr[[i]][queryHits(olaps)])
      } else {
        EWSR1_fusions[[i]] <- vcf_gr[[i]][queryHits(olaps)]
      }
    } else {
    	if (i==1) {
        EWSR1_fusions <- list(NA)
      } else {
        EWSR1_fusions[[i]] <- NA
      }
    }

  } else {
  	if (i==1) {
      EWSR1_fusions <- list(NA)
    } else {
      EWSR1_fusions[[i]] <- NA
    }
  }

  # remove duplicates:
  EWSR1_fusions[[i]] <- EWSR1_fusions[[i]][!duplicated(start(EWSR1_fusions[[i]]))]

  # remove chr22 to chr22 fusions:
  EWSR1_fusions[[i]] <- EWSR1_fusions[[i]][grep("chr22", EWSR1_fusions[[i]]$join, invert = TRUE)]

}
names(EWSR1_fusions) <- samplenames

# remove NAs:
EWSR1_fusions <- EWSR1_fusions[!is.na(EWSR1_fusions)]

for (no in c("016", "018", "031", "014", "021", "005", "019", "025", "007", "006", "012", "013", "008", "020", "015")) {
  print(EWSR1_fusions[grep(no, names(EWSR1_fusions))])
}

# count ranges:
sapply(EWSR1_fusions, length)

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
  
    if (i==1) {
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
    join_ranges[[i]] <- NA
  }  

}
names(join_ranges) <- names(EWSR1_fusions)

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

GOI_fusions <- lapply(GOI, function(x) {

  for (i in 1:length(join_ranges)) {
  
    print(i)
  
    if (i==1) {
      fusions <- list(find_fusions(join_ranges[[i]], x))
    } else {
      fusions[[i]] <- find_fusions(join_ranges[[i]], x)
    }
  
  }
  names(fusions) <- names(join_ranges)

  return(fusions)

})

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

######

# define FLI1 co-ords:
FLI1 <- GRanges(
  seqnames = "chr11",
  ranges = IRanges(start = 128554430, end = 128685162),
  strand = "*"
)


for (i in 1:length(join_ranges)) {

  print(i)

  if (i==1) {
    EWSR1_FLI1_fusions <- list(find_fusions(join_ranges[[i]], FLI1))
  } else {
    EWSR1_FLI1_fusions[[i]] <- find_fusions(join_ranges[[i]], FLI1)
  }

}
names(EWSR1_FLI1_fusions) <- names(EWSR1_fusions)

# define ETV1 co-ords:
ETV1 <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(start = 13928854, end = 14033050),
  strand = "*"
)

for (i in 1:length(join_ranges)) {

  print(i)

  if (i==1) {
    EWSR1_ETV1_fusions <- list(find_fusions(join_ranges[[i]], ETV1))
  } else {
    EWSR1_ETV1_fusions[[i]] <- find_fusions(join_ranges[[i]], ETV1)
  }

}
names(EWSR1_ETV1_fusions) <- names(EWSR1_fusions)

# define ERG co-ords:
ERG <- GRanges(
  seqnames = "chr21",
  ranges = IRanges(start = 39737183, end = 40035704),
  strand = "*"
)

for (i in 1:length(join_ranges)) {

  print(i)

  if (i==1) {
    EWSR1_ERG_fusions <- list(find_fusions(join_ranges[[i]], ERG))
  } else {
    EWSR1_ERG_fusions[[i]] <- find_fusions(join_ranges[[i]], ERG)
  }

}
names(EWSR1_ERG_fusions) <- names(EWSR1_fusions)


test_ranges <- join_ranges[!is.na(join_ranges)]

chr7 <- lapply(test_ranges, function(x) x[seqnames(x) == "chr7"])
chr21 <- lapply(test_ranges, function(x) x[seqnames(x) == "chr21"])
chr17 <- lapply(test_ranges, function(x) x[seqnames(x) == "chr17"])
chr2 <- lapply(test_ranges, function(x) x[seqnames(x) == "chr2"])

