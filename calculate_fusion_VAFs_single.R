
samplename <- "409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001"
venn_cols <- c("#7C1BE2", "#1B9E77", "#EFC000FF", "blue")
min_overlap <- 19
disc_read_window <- 200

home_dir <- "/Users/torpor/clusterHome/"
#home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
func_dir <- paste0(project_dir, "scripts/functions/")
Robject_dir <- paste0(project_dir, "Rdata/", samplename, "/")
col_dir <- paste0(home_dir, "R/colour_palettes/")

system(paste0("mkdir -p ", Robject_dir))

fusion_dir <- paste0(project_dir, "results/fusions/")
bam_path <- paste0(project_dir, "results/BWA_and_picard/bams/")


####################################################################################
### 0. Load packages and functions ###
####################################################################################

library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)
library(reshape)
library(ggplot2)

create_venn <- dget(paste0(func_dir, "create_venn.R"))
fetch_reads <- dget(paste0(func_dir, "fetch_reads.R"))
find_overlapping_reads <- dget(paste0(func_dir, "find_overlapping_reads.R"))

filter_overlaps <- function (reads, min_overlap) {
  
  split_reads <- split(reads, seqnames(reads))
  split_reads <- split_reads[c("chr11", "chr22")]
  split_reads <- lapply(split_reads, function(y) {
    if (length(y) > 0) {
      if (unique(seqnames(y)) == "chr11") {
        y$start_to_fusion <- y$fusion_coord - start(y)
        y$fusion_to_end <- end(y) - y$fusion_coord
      }
      else if (unique(seqnames(y)) == "chr22") {
        y$start_to_fusion <- y$fusion_coord - start(y)
        y$fusion_to_end <- end(y) - y$fusion_coord
      }
    }
    return(y)
  })
  prefilt_reads <- c(split_reads[[1]], split_reads[[2]])
  
  return(
    prefilt_reads[
      prefilt_reads$start_to_fusion >= min_overlap & 
        prefilt_reads$fusion_to_end >= min_overlap
    ]
  )
  
}

fetch_mate_gap <- dget(paste0(func_dir, "fetch_mate_gap.R"))
find_spanning_discordant <- dget(paste0(func_dir, "find_spanning_discordant.R"))

plot_cols <- read.table(
  paste0(col_dir, "labelled_colour_palette.txt"),
  sep = "\t",
  header = F,
  comment.char = "",
  fill = TRUE
)$V1
plot_cols <- plot_cols[c(1:3, 5, 4, 6:length(plot_cols))]


####################################################################################
### 1. Load and filter data ###
####################################################################################

# load in fusions:
fusions <- readRDS(paste0(fusion_dir, "EWSR1_GOI_fusions.Rdata"))
fusions <- fusions$collapsed$high_conf_bp$GOI_fusions$FLI1
fusions <- fusions[["409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001"]]

# split into chr22 and chr11 fusion positions:
chr22_fusions <- GRanges(
  seqnames = fusions$join_chr,
  ranges = IRanges(start = fusions$join_coord, end = fusions$join_coord),
  strand = "*",
  join_chr = seqnames(fusions),
  join_coord = start(fusions)
)
chr11_fusions <- fusions

# define file types to import:
bamtypes <- list(
  concordant_pairs = "consensus.concordant.pairs.bam", 
  discordant_pairs = "consensus.discordant.bam",
  split_supp = "consensus.split.bam"
)

# set up scanBamParam to filter out unmapped reads:
param <- ScanBamParam(
  flag = scanBamFlag(isUnmappedQuery = F),
  what = c(
    "rname", "pos", "qwidth", "strand", 
    "qname", "flag", "mapq", "cigar",
    "seq", "qual"
  )
)

if (!file.exists(paste0(Robject_dir, "VAF_calculation_reads.Rdata"))) {
  
  # load in bam:
  unfilt_bam <- lapply(bamtypes, function(x) {
    
    bam_obj <- scanBam(
      paste0(bam_path, samplename, "/", samplename, ".", x),
      param = param
    )
    
    # convert to GRanges:
    gr <- GRanges(
      seqnames = bam_obj[[1]]$rname,
      ranges = IRanges(
        start = bam_obj[[1]]$pos, 
        width = bam_obj[[1]]$qwidth
      ),
      strand = bam_obj[[1]]$strand,
      qname = bam_obj[[1]]$qname,
      flag = bam_obj[[1]]$flag,
      mapq = bam_obj[[1]]$mapq,
      cigar = bam_obj[[1]]$cigar,
      seq = bam_obj[[1]]$seq,
      qual= bam_obj[[1]]$qual
    )
    
    return(gr)
    
  })
  
  # calculate read numbers:
  read_numbers <- data.frame(
    unfiltered = sapply(unfilt_bam, length)
  )
  read_numbers["split_pairs",] = NA
  read_numbers["non_split_concordant_pairs",] = NA

  # filter all bams:
  # initiate cluster:
  cl <- makeCluster(3)
  clusterExport(
    cl, varlist = c("unfilt_bam")
  )

  temp_bam <- parLapply(cl, unfilt_bam, function(x) {

  	# remove multimapping reads (mapq score = 0):
    mmappers <- unique(x$qname[x$mapq == 0])
    res <- x[!(x$qname %in% mmappers)]
    
    # calculate read numbers:
    mmappers_removed = length(res)
  
    return(
      list(
        bam = res,
        read_no = mmappers_removed
      )
    )

  })

  stopCluster(cl)
  
  # update read number record:
  read_numbers$mmappers_removed <- c(
    sapply(temp_bam, function(x) {
      return(x$read_no)
    }),
    split_pairs = NA,
    non_split_concordant_pairs = NA
  )
  
  temp_bam <- lapply(temp_bam, function(x) return(x$bam))

  # filter concordant and discordant pairs:
  paired_bam <- temp_bam[names(temp_bam) %in% c("concordant_pairs", "discordant_pairs")]

  # initiate cluster:
  cl <- makeCluster(3)
  clusterExport(
    cl, varlist = c("paired_bam", "func_dir")
  )

  filt_bam <- lapply(paired_bam, function(x) {
    
    # remove reads with >1 supplementary alignment:
    spl <- split(x, x$qname)
    filt_spl <- spl[
      sapply(spl, function(y) {
        return(length(which(y$flag >= 2000)) <= 1)
      })
    ]
    x <- unlist(filt_spl)
    
    # record read numbers:
    read_no <- list(too_many_supp_removed = length(x))

  	# remove reads with only supplementary alignments:
    spl <- split(x, x$qname)

    filt_spl <- spl[
      !(
        sapply(spl, function(y) {
          return(all(y$flag >= 2000))
        })
      )
    ]
    res <- unlist(filt_spl)
    
    # record read numbers:
    read_no$only_supp_removed <- length(res)

    # remove unpaired reads:
    fetch_reads <- dget(paste0(func_dir, "fetch_reads.R"))
  	singles_vs_pairs <- fetch_reads(res)
  	
  	# record read numbers:
  	read_no$unpaired_removed = length(singles_vs_pairs$pairs)

    return(
      list(
        bam = singles_vs_pairs$pairs,
        read_no = unlist(read_no)
      )
    )

  })
  
  stopCluster(cl)
  
  # add read numbers to record:
  temp_read_no <- rbind(
    as.data.frame(
      t(
        sapply(filt_bam, function(x) {
          return(x$read_no)
        })
      )
    ),
    data.frame(
      too_many_supp_removed = c(
        read_numbers$mmappers_removed[rownames(read_numbers) == "split_supp"], 
        NA, NA
      ), 
      only_supp_removed = c(
        read_numbers$mmappers_removed[rownames(read_numbers) == "split_supp"], 
        NA, NA
      ),
      unpaired_removed = c(
        read_numbers$mmappers_removed[rownames(read_numbers) == "split_supp"], 
        NA, NA
      )
    )
  )
  rownames(temp_read_no)[3:5] <- c(
    "split_supp", "split_pairs", "non_split_concordant_pairs"
  )
  
  # combine with other read counts:
  read_numbers <- cbind(read_numbers, temp_read_no)

  # combine read objects into list:
  bam <- list(
  	concordant_pairs = filt_bam$concordant_pairs$bam,
  	discordant_pairs = filt_bam$discordant_pairs$bam,
  	split_supp = temp_bam$split_supp
  )

  # remove split supps not in either concordant or discordant pair elements:
  bam$split_supp <- bam$split_supp[
    bam$split_supp$qname %in% c(
      bam$concordant_pairs$qname,
      bam$discordant_pairs$qname
  	)
  ]
  
  # fetch split primary pairs:
  bam$split_pairs <- c(
  	bam$concordant_pairs[bam$concordant_pairs$qname %in% bam$split_supp$qname],
  	bam$discordant_pairs[bam$discordant_pairs$qname %in% bam$split_supp$qname]
  )
  
  # record read numbers:
  read_numbers$unmatched_split_removed <- c(
    sapply(bam, length),
    non_split_concordant_pairs = NA
  )

  # remove split pairs from discordant pairs:
  bam$discordant_pairs <- bam$discordant_pairs[
  	!(bam$discordant_pairs$qname %in% bam$split_pairs$qname)
  ]
  
  # record read numbers:
  read_numbers$split_removed_from_discordant <- c(
    sapply(bam, length),
    non_split_concordant_pairs = NA
  )
  
  # fetch non_split_concordant_pairs:
  bam$non_split_concordant_pairs <- bam$concordant_pairs[
    !(bam$concordant_pairs$qname %in% bam$split_pairs$qname)
  ]
  
  # record read numbers:
  read_numbers$split_removed_from_discordant <- sapply(bam, length)
  
  # remove pairs not mapped to either chr11 or 22:
  chr_filt_bam <- c(
    list(
      concordant_pairs = bam$concordant_pairs,
      discordant_pairs = bam$discordant_pairs,
      non_split_concordant_pairs = bam$non_split_concordant_pairs
    ),
    list(
      all_split = c(
        bam$split_pairs,
        bam$split_supp
      )
    )
  )
  spl <- lapply(chr_filt_bam, function(x) split(x, x$qname))
  
  # initiate cluster:
  cl <- makeCluster(7)
  clusterExport(
    cl, varlist = c("spl")
  )
  
  # remove concordant pairs not mapped to either chr11 or 22:
  spl$concordant_pairs <- spl$concordant_pairs[
    unlist(
      parLapply(cl, spl$concordant_pairs, function(x) {
        all(seqnames(x) %in% c("chr11", "chr22"))
      })
    )
  ]

  # remove non split concordant pairs not mapped to chr11 and 22:
  spl$non_split_concordant_pairs <- spl$non_split_concordant_pairs[
    unlist(
      parLapply(cl, spl$non_split_concordant_pairs, function(x) {
        all(seqnames(x) %in% c("chr11", "chr22"))
      })
    )
  ]
  
  # remove discordant pairs not mapped to chr11 and 22:
  spl$discordant_pairs <- spl$discordant_pairs[
    unlist(
      parLapply(cl, spl$discordant_pairs, function(x) {
        "chr11" %in% seqnames(x) & "chr22" %in% seqnames(x)
      })
    )
  ]
  
  # remove split reads not mapped to chr11 and 22:
  spl$all_split <- spl$all_split[
    unlist(
      parLapply(cl, spl$all_split, function(x) {
        "chr11" %in% seqnames(x) & "chr22" %in% seqnames(x)
      })
    )
  ]
  
  stopCluster(cl)
  
  # merge each set of reads and separate split supp and primary:
  fusion_chr_bam <- lapply(spl, function(x) unlist(x))
  fusion_chr_bam$split_supp <- bam$split_supp[
    bam$split_supp$qname %in% fusion_chr_bam$all_split$qname
  ]
  fusion_chr_bam$split_pairs <- bam$split_pairs[
    bam$split_pairs$qname %in% fusion_chr_bam$all_split$qname
  ]
  fusion_chr_bam <- fusion_chr_bam[
    !(names(fusion_chr_bam) %in% "all_split")
  ]
  
  # remove rownames:
  fusion_chr_bam <- lapply(fusion_chr_bam, function(x) {
    names(x) <- NULL
    return(x)
  })

  # record read numbers:
  temp_no <- sapply(fusion_chr_bam, length)
  read_numbers$fusion_chr_only <- temp_no[
    match(rownames(read_numbers), names(temp_no))
  ]

  saveRDS(unfilt_bam, paste0(Robject_dir, "unfiltered_VAF_calculation_reads.Rdata"))
  saveRDS(fusion_chr_bam, paste0(Robject_dir, "VAF_calculation_reads.Rdata"))
  saveRDS(read_numbers, paste0(Robject_dir, "VAF_calculation_read_nos.Rdata"))
  
} else {

  unfilt_bam <- readRDS(paste0(Robject_dir, "unfiltered_VAF_calculation_reads.Rdata"))
  fusion_chr_bam <- readRDS(paste0(Robject_dir, "VAF_calculation_reads.Rdata"))
  read_numbers <- readRDS(paste0(Robject_dir, "VAF_calculation_read_nos.Rdata"))

}

# create venn diagram of filtered vs unfiltered reads:
filt_vs_unfilt <- list(
  unfiltered = c(unfilt_bam$concordant_pairs, unfilt_bam$discordant_pairs), 
  filtered = c(fusion_chr_bam$concordant_pairs, fusion_chr_bam$discordant_pairs)
)
filt_vs_unfilt_venn <- create_venn(filt_vs_unfilt, venn_cols)

# create venn diagram of split concordant pairs:
concordant_split <- fusion_chr_bam[
  names(fusion_chr_bam) %in% c(
    "concordant_pairs", "split_pairs", "split_supp"
  )
]
concordant_split_venn <- create_venn(concordant_split, venn_cols)

# create venn diagram of non-split concordant pairs:
concordant_non_split <- fusion_chr_bam[
  names(fusion_chr_bam) %in% c(
    "concordant_pairs", "non_split_concordant_pairs"
  )
]
concordant_non_split_venn <- create_venn(concordant_non_split, venn_cols)

# create venn diagram of discordant pairs:
discordant <- fusion_chr_bam[
  names(fusion_chr_bam) %in% c(
    "discordant_pairs", "split_pairs", "split_supp"
  )
]
discordant_venn <- create_venn(discordant, venn_cols)

# create barplot of read filtering:
read_numbers$type <- gsub("_", " ", rownames(read_numbers))
colnames(read_numbers) <- gsub("_", " ", colnames(read_numbers))
plot_df <- melt(read_numbers)
plot_df$type <- factor(
  plot_df$type, 
  levels = c(
    "concordant pairs", "discordant pairs", "split supp", 
    "split pairs", "non split concordant pairs"
  )
)

# plot:
p <- ggplot(plot_df, aes(x=variable, y=value, fill=type))
p <- p + geom_bar(stat="identity", position = "dodge")
p <- p + scale_y_continuous(trans='log10')
p <- p + scale_fill_manual(values=plot_cols)
p <- p + ylab("No. reads (log10)")
p <- p + theme_cowplot(12)
p <- p + theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x=element_blank()
)


####################################################################################
### 2. find breakpoint-overlapping reads: ###
####################################################################################

# create empty list to be filled:
fusion_reads <- list(
  spanning = NULL,
  overlapping = NULL
)

# find overlaps of non-split concordant pairs with chr22 breakpoint coords, 
# add to non_supporting_reads$overlapping:
non_supporting_reads <- fusion_reads

# split concordant reads by qname:
spl <- split(
  fusion_chr_bam$non_split_concordant_pairs, 
  fusion_chr_bam$non_split_concordant_pairs$qname
)

# initiate cluster:
cl <- makeCluster(7)
clusterExport(
  cl, varlist = c("spl")
)

# find overlaps with fusion breakpoint:
overlapping_non_supporting_reads <- parLapply(
  cl,
  spl,
  find_overlapping_reads,
  fusions = chr22_fusions,
  chromosome = "chr22"
)

stopCluster(cl)

non_supporting_reads$overlapping <- unlist(
  as(
    overlapping_non_supporting_reads[
      sapply(overlapping_non_supporting_reads, function(x) !is.null(x))
    ],
    "GRangesList"
  )
)
names(non_supporting_reads$overlapping) <- NULL

# find overlaps of split primary reads with chr22 breakpoint coords, 
# add to supporting_reads$overlapping:
supporting_reads <- fusion_reads

spl <- split(fusion_chr_bam$split_pairs, fusion_chr_bam$split_pairs$qname)

# initiate cluster:
cl <- makeCluster(7)
clusterExport(
  cl, varlist = c("spl")
)

# find overlaps with fusion breakpoint:
overlapping_supporting_reads <- parLapply(
  cl,
  spl,
  find_overlapping_reads,
  fusions = chr22_fusions,
  chromosome = "chr22"
)

stopCluster(cl)

supporting_reads$overlapping <- unlist(
  as(
    overlapping_supporting_reads[
      sapply(overlapping_supporting_reads, function(x) !is.null(x))
    ],
    "GRangesList"
  )
)
names(supporting_reads$overlapping) <- NULL

# find overlaps of split primary reads with chr11 breakpoint coords, 
# add to supporting_reads$overlapping:

# initiate cluster:
cl <- makeCluster(7)
clusterExport(
  cl, varlist = c("spl")
)

# find overlaps with fusion breakpoint:
overlapping_supporting_reads <- parLapply(
  cl,
  spl,
  find_overlapping_reads,
  fusions = chr11_fusions,
  chromosome = "chr11"
)

stopCluster(cl)

supporting_reads$overlapping <- c(
  supporting_reads$overlapping,
  unlist(
    as(
      overlapping_supporting_reads[
        sapply(overlapping_supporting_reads, function(x) !is.null(x))
      ],
      "GRangesList"
    )
  )
)
names(supporting_reads$overlapping) <- NULL

# deduplicate supporting reads:
spl <- split(supporting_reads$overlapping, supporting_reads$overlapping$qname)
overlapping_supporting_reads <- lapply(spl, function(x) {
  return(x[!duplicated(x)])
})

supporting_reads$overlapping <- unlist(
  as(
    overlapping_supporting_reads[sapply(overlapping_supporting_reads, function(x) !is.null(x))],
    "GRangesList"
  )
)
names(supporting_reads$overlapping) <- NULL

# deduplicate non_supporting reads:
spl <- split(non_supporting_reads$overlapping, non_supporting_reads$overlapping$qname)
overlapping_non_supporting_reads <- lapply(spl, function(x) {
  return(x[!duplicated(x)])
})

non_supporting_reads$overlapping <- unlist(
  as(
    overlapping_non_supporting_reads[sapply(overlapping_non_supporting_reads, function(x) !is.null(x))],
    "GRangesList"
  )
)
names(non_supporting_reads$overlapping) <- NULL

# filter overlapping reads without at least min_overlap on both sides of fusion:
non_supporting_reads$overlapping <- filter_overlaps(
  reads = non_supporting_reads$overlapping, 
  min_overlap = min_overlap
)

supporting_reads$overlapping <- filter_overlaps(
  reads = supporting_reads$overlapping, 
  min_overlap = min_overlap
)


####################################################################################
### 3. Find breakpoint-spanning reads ###
####################################################################################

# find breakpoint-spanning reads from non-split non-discordant reads, 
# add to spanning$non_supporting:

# fetch read gaps (read + gap coords):
# split by qname:
spl <- split(
  fusion_chr_bam$non_split_concordant_pairs, 
  fusion_chr_bam$non_split_concordant_pairs$qname
)

# initiate cluster:
cl <- makeCluster(7)
clusterExport(
  cl, varlist = c("spl")
)

# fetch read gaps:
non_split_concordant_gaps <- parLapply(
  cl, spl, fetch_mate_gap
)

stopCluster(cl)

# merge into granges:
non_split_concordant_gaps <- unlist(
  as(non_split_concordant_gaps, "GRangesList")
)
names(non_split_concordant_gaps) <- NULL

# find overlaps of gaps with fusions:
spanning_non_supporting_gaps <- find_overlapping_reads(
  fusions = chr22_fusions, 
  reads = non_split_concordant_gaps,
  chromosome = "chr22"
)

# fetch corresponding read pairs:
non_supporting_reads$spanning <- fusion_chr_bam$non_split_concordant_pairs[
  fusion_chr_bam$non_split_concordant_pairs$qname %in% 
    spanning_non_supporting_gaps$qname
]
m <- match(
  non_supporting_reads$spanning$qname, spanning_non_supporting_gaps$qname
)
non_supporting_reads$spanning$chr22_fusion_coord <- 
  spanning_non_supporting_gaps$chr22_fusion_coord[m]

# find breakpoint-spanning reads from discordant reads, 
# add to spanning$supporting:
spl <- split(
  fusion_chr_bam$discordant_pairs, 
  fusion_chr_bam$discordant_pairs$qname
)
spanning_supporting_reads <- lapply(
  spl, find_spanning_discordant, 
  chr11_fusions, chr22_fusions, disc_read_window
)

# merge to granges object:
supporting_reads$spanning <- unlist(
  as(
    spanning_supporting_reads[
      sapply(spanning_supporting_reads, function(x) !is.null(x))
    ], "GRangesList"
  )
)
names(supporting_reads$spanning) <- NULL


####################################################################################
### 4. Calculate proportions of non-supporting vs supporting read pairs for each
# sample ###
####################################################################################

# check all reads are in pairs:
supporting_check <- sapply(supporting_reads, function(x) {
  spl <- split(x, x$qname)
  return(
    all(
      sapply(spl, function(y) {
        length(y) == 2
      })
    )
  )
})
print(paste0("Are all supporting reads in pairs? ", all(supporting_check)))

non_supporting_check <- sapply(non_supporting_reads, function(x) {
  spl <- split(x, x$qname)
  return(
    all(
      sapply(spl, function(y) {
        length(y) == 2
      })
    )
  )
})
print(
  paste0("Are all non-supporting reads in pairs? ", all(non_supporting_check))
)

# remove non-primary read mate from each pair in order to split by fusion. 
# pairs with one mate supporting one breakpoint, one supporting another
# will be counted towards the two breakpooints when split:
remove_pairs <- function(reads) {
  
  # split reads by name:
  spl <- split(reads, reads$qname)
  
  res <- lapply(spl, function(x) {
    
    # remove read with no overlapping fusion coord:
    filt_read <- x[!is.na(x$fusion_coord)]
    
    # if there are still 2 reads, remove the one mapped to chromosome 11:
    if (length(filt_read) > 1) {
      filt_read <- filt_read[filt_read$fusion_chr == "chr11"]
    }
    
    # if there are still 2 reads and they have the same fusion coordinate, remove 
    # the first one:
    if (length(filt_read) > 1) {
      filt_read <- filt_read[1]
    }
    
    return(filt_read)
    
  })
  
  res <- unlist(
    as(
      res, "GRangesList"
    )
  )
  names(res) <- NULL
  
  return(res)
  
}

fusion_supporting_reads <- lapply(supporting_reads, remove_pairs)



# check numbers for following - don't make sense:
fusion_non_supporting_reads <- lapply(non_supporting_reads, remove_pairs)


# remove non-split reads from split read pairs:
fusion_supporting_reads$overlapping <- supporting_reads$overlapping[
  !is.na(supporting_reads$overlapping$fusion_coord)
]

# remove non-split reads from non-split overlapping read pairs:
fusion_non_supporting_reads$overlapping <- non_supporting_reads$overlapping[
  !is.na(non_supporting_reads$overlapping$fusion_coord)
]

# remove chr11 reads from discordant read pairs:
fusion_supporting_reads$spanning <- supporting_reads$spanning[
  supporting_reads$spanning$fusion_chr != "chr22"
]

# remove chr11 reads from concordant spanning read pairs:
fusion_supporting_reads$spanning <- supporting_reads$spanning[
  supporting_reads$spanning$fusion_chr != "chr22"
]


######
all_supporting <- unlist(as(fusion_supporting_reads, "GRangesList"))
######

# split non_supporting read mates into fusions:
fusion_non_supporting_reads <- lapply(
  non_supporting_reads, function(x) split(x, x$fusion_coord)
)

  
  res <- unlist(
    as(
      combined[
        sapply(combined, function(x) !is.null(x))
      ], "GRangesList"
    )
  )
  names(res) <- NULL
  
  return(res)
  
})

# keep breakpoint with greatest number of supporting reads:
lapply(supporting_reads, function(x) split(x, x$chr22_fusion_coord))

# add overlapping and spanning reads for supporting and non-supporting filtered
# reads:

# calculate supporting vs non-supporting proportion:

