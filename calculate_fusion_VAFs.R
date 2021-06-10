
filter_multimappers <- TRUE
filter_max_supp_align <- 1  # make NA to switch off filter
filter_singles <- TRUE
filter_supp_only <- TRUE
min_overlap <- 19
venn_cols <- c("#7C1BE2", "#1B9E77", "#EFC000FF", "blue", "brown")

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
#home_dir <- "/Users/jamestorpy/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
func_dir <- paste0(project_dir, "scripts/functions/")
Robject_dir <- paste0(project_dir, "Rdata/")

system(paste0("mkdir -p ", Robject_dir))

fusion_dir <- paste0(project_dir, "results/fusions/")
bam_path <- paste0(project_dir, "results/BWA_and_picard/bams/")

####################################################################################
### 0. Load packages and functions ###
####################################################################################

library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)
library(ggvenn)

find_overlapping_reads <- dget(paste0(func_dir, "find_overlapping_reads.R"))
filter_overlaps <- dget(paste0(func_dir, "filter_overlaps.R"))
create_venn <- dget(paste0(func_dir, "create_venn.R"))
fetch_reads <- dget(paste0(func_dir, "fetch_reads.R"))


####################################################################################
### 1. Load and filter data ###
####################################################################################

# load in fusions:
fusions <- readRDS(paste0(fusion_dir, "EWSR1_GOI_fusions.Rdata"))
fusions <- fusions$collapsed$high_conf_bp$GOI_fusions$FLI1

# define samplenames:
samplenames <- names(fusions)

# split into chr22 and chr11 fusion positions:
chr22_fusions <- lapply(fusions, function(x) {
  return(
    GRanges(
      seqnames = x$join_chr,
      ranges = IRanges(start = x$join_coord, end = x$join_coord),
      strand = "*",
      join_chr = seqnames(x),
      join_coord = start(x)
    )
  )
})
chr11_fusions <- fusions

# define file types to import:
bamtypes <- list(
  all = "consensus.bam", 
  split = "consensus.split.bam", 
  discordant = "consensus.discordant.bam"
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

if (!file.exists(paste0(Robject_dir, "VAF_calculation_bams.Rdata"))) {
  
  # load in bam, split bam and discordant bam as GRanges:
  unfilt_bams <- lapply(samplenames, function(x) {
    
    # import bams:
    bam_temp <- lapply(bamtypes, function(y) {
      
      bam_obj <- scanBam(
        paste0(bam_path, x, "/", x, ".", y),
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
      
      ######
      #gr <- gr[1:floor(length(gr)/10)]
      ######
      
      return(gr)

    })
    
  })
  names(unfilt_bams) <- samplenames
  
  # filter bams:
  bams <- lapply(unfilt_bams, function(x) {
    
    return(
      lapply(x, function(y) {
        
        if (filter_multimappers) {
          # remove multimapping reads (mapq score = 0):
          mmappers <- unique(y$qname[y$mapq == 0])
          y <- y[!(y$qname %in% mmappers)]
        }
        
        if (!is.na(filter_max_supp_align)) {
          # remove reads with >1 supplementary alignment:
          spl <- split(y, y$qname)
          filt_spl <- spl[
            sapply(spl, function(z) {
              return(length(which(z$flag >= 2000)) <= filter_max_supp_align)
            })
          ]
        }
        
        if (filter_supp_only) {
          # remove reads with only supplementary alignments:
          spl <- split(y, y$qname)
          filt_spl <- spl[
            sapply(spl, function(z) {
              return(length(which(z$flag < 2000)) <= 1)
            })
          ]
        }
        
        return(unlist(filt_spl))
        
      })
    )
    
  })
  saveRDS(bams, paste0(Robject_dir, "filtered_bams.Rdata"))
  
  ####################################################################################
  ### 2. Split reads into groups and plot ###
  ####################################################################################
  
  # fetch discarded reads for venn below:
  for (i in 1:length(unfilt_bams)) {
    bams[[i]]$all_discarded <- unfilt_bams[[i]]$all[
      !(unfilt_bams[[i]]$all$qname %in% bams[[i]]$all$qname)
    ]
  }
  
  for (i in 1:length(bams)) {
    
    writeLines("\n")
    
    # remove discordant reads from split reads and vice versa:
    print(
      paste0("Removing discordant from split reads for ", names(bams)[i], "...")
    )
    
    bams[[i]]$split <- bams[[i]]$split[
      !(bams[[i]]$split$qname %in% bams[[i]]$discordant$qname)
    ]
    bams[[i]]$discordant <- bams[[i]]$discordant[
      !(bams[[i]]$discordant$qname %in% bams[[i]]$split$qname)
    ]
    
    print(
      paste0("Fetching non-split non-discordant reads of ", names(bams)[i], "...")
    )
    
    # combine discordant and split qnames:
    discordant_and_split <- c(
      bams[[i]]$discordant$qname,
      bams[[i]]$split$qname
    )
    
    # fetch non discordant/split reads:
    bams[[i]]$non_discordant_or_split <- bams[[i]]$all[
      !(bams[[i]]$all$qname %in% discordant_and_split)
    ]
    
    print(
      paste0("Fetching split primary reads of ", names(bams)[i], "...")
    )
    
    # create split_primary GRanges:
    bams[[i]]$split_primary <- bams[[i]]$all[
      bams[[i]]$all$qname %in% bams[[i]]$split$qname
    ]
    bams[[i]]$split_primary <- bams[[i]]$split_primary[
      as.numeric(bams[[i]]$split_primary$flag) < 2000
    ]
    
    print(
      paste0("Separating paired and unpaired reads of ", names(bams)[i], "...")
    )
    
    # subset to separate reads:
    bams_sub <- bams[[i]][
      names(bams[[i]]) %in% c(
        "all", "split", "split_primary", "discordant", "non_discordant_or_split"
      )
    ]
    
    # isolate paired and single reads:
    bams_sub <- lapply(bams_sub, fetch_reads)
    
    # add all to central list:
    bams[[i]] <- c(
      bams[[i]],
      list(
        pairs = list(
          all = bams_sub$all$pairs,
          split = bams_sub$split$pairs,
          discordant = bams_sub$discordant$pairs,
          non_discordant_or_split = bams_sub$non_discordant_or_split$pairs,
          split_primary = bams_sub$split_primary$pairs
        )
      ),
      list(
        singles = list(
          all = bams_sub$all$singles,
          split = bams_sub$split$singles,
          discordant = bams_sub$discordant$singles,
          non_discordant_or_split = bams_sub$non_discordant_or_split$singles,
          split_primary = bams_sub$split_primary$singles
        )
      )
    )
    
  }
  
  saveRDS(bams, paste0(Robject_dir, "VAF_calculation_bams.Rdata"))
  
} else {
  bams <- readRDS(paste0(Robject_dir, "VAF_calculation_bams.Rdata"))
}

# create venn diagram of filtered vs unfiltered reads:
for (i in 1:length(bams)) {
  
  if (i==1) {
    filt_vs_unfilt <- list(
      list(
        unfiltered = unfilt_bams[[i]]$all, 
        filtered = bams[[i]]$all, 
        discarded = bams[[i]]$all_discarded
      )
    )
  } else {
    filt_vs_unfilt[[i]] <- list(
      unfiltered = unfilt_bams[[i]]$all, 
      filtered = bams[[i]]$all, 
      discarded = bams[[i]]$all_discarded
    )
  }
  
}
names(filt_vs_unfilt) <- names(bams)

filt_vs_unfilt_venns <- lapply(filt_vs_unfilt, create_venn, venn_cols)

# create venn diagram of paired vs single reads:
for (i in 1:length(bams)) {
  
  if (i==1) {
    pairs_vs_singles <- list(
      list(
        all = bams[[i]]$all, 
        pairs = bams[[i]]$pairs$all, 
        singles = bams[[i]]$singles$all
      )
    )
  } else {
    pairs_vs_singles[[i]] <- list(
      all = bams[[i]]$all, 
      pairs = bams[[i]]$pairs$all, 
      singles = bams[[i]]$singles$all
    )
  }
  
  # mark supp reads as distinct:
  pairs_vs_singles[[i]] <- lapply(pairs_vs_singles[[i]], function(x) {
    x$qname[x$flag >= 2000] <- paste0(
      x$qname[x$flag >= 2000],
      "_supp"
    )
    return(x)
  })
  
}
names(pairs_vs_singles) <- names(bams)

pairs_vs_singles_venns <- lapply(pairs_vs_singles, create_venn, venn_cols)

# create venn diagram of paired reads:
for (i in 1:length(bams)) {
  
  if (i==1) {
    all_pairs <- list(
      list(
        all = bams[[i]]$pairs$all, 
        split = bams[[i]]$pairs$split, 
        split_primary = bams[[i]]$pairs$split_primary,
        discordant = bams[[i]]$pairs$discordant,
        non_discordant_or_split = bams[[i]]$pairs$non_discordant_or_split
      )
    )
    names(all_pairs[[i]]) <- c(
      "all", "split supplementary\nalignment", "split primary\nalignment",
      "discordant", "non-discordant or split"
    )
  } else {
    all_pairs[[i]] <- list(
      all = bams[[i]]$pairs$all, 
      split = bams[[i]]$pairs$split, 
      split_primary = bams[[i]]$pairs$split_primary,
      discordant = bams[[i]]$pairs$discordant,
      non_discordant_or_split = bams[[i]]$pairs$non_discordant_or_split
    )
    names(all_pairs[[i]]) <- c(
      "all", "split supplementary\nalignment", "split primary\nalignment",
      "discordant", "non-discordant or split"
    )
  }
  
}
names(all_pairs) <- names(bams)

all_pairs_venns <- lapply(all_pairs, create_venn, venn_cols)

# create venn diagram of unpaired reads:
for (i in 1:length(bams)) {
  
  if (i==1) {
    all_singles <- list(
      list(
        all = bams[[i]]$singles$all, 
        split = bams[[i]]$singles$split, 
        split_primary = bams[[i]]$singles$split_primary,
        discordant = bams[[i]]$singles$discordant,
        non_discordant_or_split = bams[[i]]$singles$non_discordant_or_split
      )
    )
    names(all_singles[[i]]) <- c(
      "all", "split supplementary\nalignment", "split primary\nalignment",
      "discordant", "non-discordant or split"
    )
  } else {
    all_singles[[i]] <- list(
      all = bams[[i]]$singles$all, 
      split = bams[[i]]$singles$split, 
      split_primary = bams[[i]]$singles$split_primary,
      discordant = bams[[i]]$singles$discordant,
      non_discordant_or_split = bams[[i]]$singles$non_discordant_or_split
    )
    names(all_singles[[i]]) <- c(
      "all", "split supplementary\nalignment", "split primary\nalignment",
      "discordant", "non-discordant or split"
    )
  }
  
}
names(all_singles) <- names(bams)

all_singles_venns <- lapply(all_singles, create_venn, venn_cols)


####################################################################################
### 2. Isolate relevant reads and plot ###
####################################################################################

# select relevant groups only:
fbams <- lapply(bams, function(x) {
  
  if (filter_singles) {
    
    writeLines("\n")
    
    res <- list(
      split_primary = x$pairs$split_primary,
      split_supp = x$singles$split,
      discordant = x$pairs$discordant,
      non_split_or_discordant = x$pairs$non_discordant_or_split
    )
    
    # remove supplementary split reads not in primary split read pairs:
    res$split_supp <- res$split_supp[
      res$split_supp$qname %in% res$split_primary$qname
    ]
    
    print(
      paste0(
        "Does the number of supplementary split alignments (", 
        length(res$split_supp),
        ") equal half the number of primary split pair alignments (",
        length(res$split_primary), ")? ",
        length(res$split_supp) == (length(res$split_primary))/2
      )
    )
    
    return(res)
    
  } else {
    
    res <- list(
      split_primary = c(x$pairs$split_primary, x$singles$split_primary),
      split_supp = x$singles$split,
      discordant = c(x$pairs$discordant, x$singles$discordant),
      non_split_or_discordant = c(
        x$pairs$non_discordant_or_split, x$singles$non_discordant_or_split
      )
    )
    
    # remove supplementary split reads not in primary split read pairs:
    res$split_supp <- res$split_supp[
      res$split_supp$qname %in% res$split_primary$qname
    ]
    
    return(res)
    
  }
  
})

# isolate reads mapping to either chr11 or 22:
fbams <- lapply(fbams, function (x) {
  lapply(x, function(y) {
    spl <- split(y, y$qname)
    keep <- lapply(spl, function(z) {
      if (all(seqnames(z) %in% c("chr11", "chr22"))) {
        return(z)
      } else {
        return(NULL)
      }
    })
    # remove NULLs:
  })
})




######

# create venn diagrams of chr11/22 read breakdowns:
fusion_gene_venns <- lapply(fbams, create_venn, venn_cols)

chr11_22_venn_vectors <- lapply(chr11_22_read_vectors, function(x) {
  res <- x[!(names(x) %in% c("all", "non_discordant_or_split"))]
  res <- res[order(names(res))]
  names(res) <- c(
    "discordant\n(not split)", "split supplementary\nalignment", 
    "split primary\nalignment"
  )
  return(res)
})

chr11_22_read_venn <- ggvenn(
  chr11_22_venn_vectors[["409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001"]], 
  fill_color = c("#7C1BE2", "#1B9E77", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)
saveRDS(chr11_22_read_venn, paste0(Robject_dir, "chr_11_22_read_venn.Rdata"))


####################################################################################
### 2. find breakpoint-overlapping reads: ###
####################################################################################

# create empty list to be filled:
fusion_reads <- list(
  spanning = NULL,
  overlapping = NULL
)

# find overlaps of non_discordant_or_split reads with chr22 breakpoint coords, 
# add to non_supporting_reads$overlapping:
for (i in 1:length(chr22_fusions)) {
  
  if (i==1) {
    non_supporting_reads <- list(fusion_reads)
  } else {
    non_supporting_reads[[i]] <- fusion_reads
  }
  
  non_supporting_reads[[i]]$overlapping <- find_overlapping_reads(
    fusions = chr22_fusions[[i]], 
    reads = bams[[i]]$non_discordant_or_split,
    chromosome = "chr22"
  )
 
}
names(non_supporting_reads) <- names(chr22_fusions)

# find overlaps of split primary reads with chr22 breakpoint coords, 
# add to supporting_reads$overlapping:
for (i in 1:length(chr22_fusions)) {
  
  if (i==1) {
    supporting_reads <- list(fusion_reads)
  } else {
    supporting_reads[[i]] <- fusion_reads
  }
  
  supporting_reads[[i]]$overlapping <- find_overlapping_reads(
    fusions = chr22_fusions[[i]], 
    reads = bams[[i]]$split_primary,
    chromosome = "chr22"
  )
  
}
names(supporting_reads) <- names(chr22_fusions)

# find overlaps of split primary reads with chr11 breakpoint coords, 
# add to supporting_reads$overlapping:
for (i in 1:length(chr11_fusions)) {
  
  supporting_reads[[i]]$overlapping <- c(
    supporting_reads[[i]]$overlapping,
    find_overlapping_reads(
      fusions = chr11_fusions[[i]], 
      reads = bams[[i]]$split_primary,
      chromosome = "chr11"
    )
  )
  
}

# filter overlapping reads without at least min_overlap on both sides of fusion:
non_supporting_reads <- lapply(non_supporting_reads, function(x) {
  x$overlapping <- filter_overlaps(reads = x$overlapping, min_overlap = min_overlap)
  return(x)
})
supporting_reads <- lapply(supporting_reads, function(x) {
  x$overlapping <- filter_overlaps(reads = x$overlapping, min_overlap = min_overlap)
  return(x)
})


####################################################################################
### 3. Find breakpoint-spanning reads ###
####################################################################################

fetch_mate_gaps <- function(reads) {
  
  library(rtracklayer)
  library(GenomicRanges)
  
  reads <- reads$non_discordant_or_split
  
  strand(reads) <- "*"
  
  # split by qname:
  qname_split <- split(reads, reads$qname)
  
  # fetch between sequences:
  mate_gaps <- lapply(qname_split, function(z) {
    
    print(paste0("Finding gaps for read id ", z$qname, "..."))
    
    if (length(z) == 2) {
      mgap <- gaps(z)
      mgap <- mgap[start(mgap) != 1]
      return(mgap)
    } else if (length(z) == 1) {
      print(
        paste0(
          z$qname, " does not have a corresponding mate, removing..."
        )
      )
      return(NA)
    } else if (length(z) > 2) {
      print(
        paste0(z$qname, " has more than one corresponding mate, removing...")
      )
    }
    
  })
  
}

if (paste0(Robject_dir, "non_discordant_or_split_read_gaps.Rdata")) {
  
  # for each non-discordant non-split read pair, fetch sequence in between mates:
  library(parallel)
  
  # calculate the number of cores:
  no_cores <- detectCores() - 1
  
  # initialise cluster:
  cl <- makeCluster(no_cores)
  
  # pass variables to cluster:
  clusterExport(cl, varlist = c("bams", "fetch_mate_gaps"))
  
  # fetch gaps between read pairs:
  # took 4.5 hours on cluster with 5 cores:
  #system.time(non_discordant_or_split_gaps <- lapply(bams, fetch_mate_gaps))
  # took 2 hours on Rstudio with 3 cores, 1 hour 20 on cluster with 5 cores:
  system.time(non_discordant_or_split_gaps <- parLapply(cl, bams, fetch_mate_gaps)) 

  
  saveRDS(
    non_discordant_or_split_gaps, 
    paste0(Robject_dir, "non_discordant_or_split_read_gaps.Rdata")
  )
  
} else {
  
  non_discordant_or_split_gaps <- readRDS(
    paste0(Robject_dir, "non_discordant_or_split_read_gaps.Rdata")
  )
  
}

# find breakpoint-spanning reads from non-split non-discordant reads, 
# add to spanning$non_supporting:

# find breakpoint-spanning reads from discordant reads, 
# add to spanning$supporting:

# calculate distance between reads and breakpoint

# keep read pairs with one within 150 bp of chr22 breakpoint, one within 150 bp 
# of chr11 breakpoint:


####################################################################################
### 4. Calculate proportions of non-supporting vs supporting read pairs for each
# sample ###
####################################################################################

# keep breakpoint with greatest number of supporting reads:

# add overlapping and spanning reads for supporting and non-supporting filtered
# reads:

# calculate supporting vs non-supporting proportion:

