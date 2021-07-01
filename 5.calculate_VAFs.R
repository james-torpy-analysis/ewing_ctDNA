
samplename <- "409_002_D9YW9_GGACTCCT-CTCTCTAT_L001"
venn_cols <- c("#7C1BE2", "#1B9E77", "#EFC000FF", "blue")
min_overlap <- 19
disc_read_window <- 200
same_fusion_window <- 10

home_dir <- "/Users/torpor/clusterHome/"
#home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
func_dir <- paste0(project_dir, "scripts/functions/")
result_dir <- paste0(project_dir, "results/")
Robject_dir <- paste0(result_dir, "VAF_calculation/", samplename, "/Rdata/")
table_dir <- paste0(result_dir, "VAF_calculation/", samplename, "/tables/")
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", table_dir))

fusion_dir <- paste0(result_dir, "fusions/", samplename, "/")
bam_path <- paste0(result_dir, "BWA_and_picard/bams/")
col_dir <- paste0(home_dir, "R/colour_palettes/")


####################################################################################
### 0. Load packages and functions ###
####################################################################################

library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)
library(reshape)
library(ggplot2)
library(cowplot)
library(ggvenn)

create_venn <- dget(paste0(func_dir, "create_venn.R"))
find_overlapping_reads <- dget(paste0(func_dir, "find_overlapping_reads.R"))

filter_overlaps <- function(reads, min_overlap) {
  
  # split by chromosome:
  split_reads <- split(reads, seqnames(reads))
  split_reads <- split_reads[c("chr11", "chr22")]
  
  # calculate lengths from start of read to fusion, and from fusion to end:
  split_reads <- lapply(split_reads, function(y) {
    if (length(y) > 0) {
      if (unique(seqnames(y)) == "chr11") {
        y$start_to_fusion <- y$fusion_coord - start(y)
        y$fusion_to_end <- end(y) - y$fusion_coord
      } else if (unique(seqnames(y)) == "chr22") {
        y$start_to_fusion <- y$fusion_coord - start(y)
        y$fusion_to_end <- end(y) - y$fusion_coord
      }
    }
    return(y)
  })
  
  # concatentate reads back together:
  prefilt_reads <- c(split_reads[[1]], split_reads[[2]])
  
  # split by qname:
  prefilt_split <- split(prefilt_reads, prefilt_reads$qname)
  
  # filter out reads without overlaps of at least min_overlap on both sides of
  # fusion (need at least 1/2 read to satisfy this condition):
  filt_split <- lapply(prefilt_split, function(x) {
    condition_vec <- !(
      x$start_to_fusion >= min_overlap & x$fusion_to_end >= min_overlap
    )
    condition_vec[is.na(condition_vec)] <- FALSE
    if (all(condition_vec)) {
      return(NULL)
    } else {
      return(x)
    }
  })
  
  res <- unlist(
    as(
      filt_split[sapply(filt_split, function(x) !is.null(x))],
      "GRangesList"
    )
  )
  names(res) <- NULL
  
  return(res)
  
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
### 1. Load data ###
####################################################################################

# load in fusions:
fusions <- readRDS(paste0(fusion_dir, "EWSR1_GOI_fusions.Rdata"))

fusions <- fusions$high_conf$true_positives$fusions$FLI1
# make chr22 coord main fusion coord:
fusions <- lapply(fusions, function(x) {
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

# combine fusions with <= same_fusion_window bp difference in coord:
fusions <- lapply(fusions, function(x) {
  
  x$remove <- FALSE
  
  for (i in 1:length(x)) {
    
    if (!(x$remove[i])) {
      
      chr22_window <- (start(x[i])-same_fusion_window):(start(x[i])+same_fusion_window)
      chr11_window <- (x[i]$join_coord-same_fusion_window):(x[i]$join_coord+same_fusion_window)
      
      # check each other fusion to see whether both chr11 and 22 coord is close
      # enough to those of another fusion to justify removing one of them:
      for (j in 1:length(x)) {
        
        if (j!=i & start(x)[j] %in% chr22_window & 
            x$join_coord[j] %in% chr11_window) {
          x$remove[j] <- TRUE
        }
        
      }
      
    }
    
  }
    
  # remove those marked as the same as another fusion:
  res <- x[!(x$remove)]
  mcols(res) <- subset(res, select = -remove)
  
  return(res)
  
})

# load filtered reads:
filtered_reads <- readRDS(paste0(Robject_dir, "VAF_calculation_reads.Rdata"))


####################################################################################
### 2. find breakpoint-overlapping reads: ###
####################################################################################

if (!file.exists(paste0(Robject_dir, "high_conf_fusion_reads.Rdata"))) {
  
  fusions <- fusions$high_conf
    
  # create empty list to be filled:
  temp_list <- list(
    spanning = NULL,
    overlapping = NULL
  )
  
  for (i in 1:length(fusions)) {
    if (i==1) {
      all_reads <- list(
        list(
          non_supporting = temp_list,
          supporting = temp_list
        )
      )
    } else {
      all_reads[[i]] <- list(
        non_supporting = temp_list,
        supporting = temp_list
      )
    }
  }
  names(all_reads) <- start(fusions)
  
  # find overlaps of non-split concordant pairs with chr22 breakpoint coords, 
  # add to non_supporting$overlapping:
  
  # split concordant reads by qname:
  spl <- split(
    filtered_reads$non_split_concordant_pairs, 
    filtered_reads$non_split_concordant_pairs$qname
  )
  
  # initiate cluster:
  cl <- makeCluster(7)
  clusterExport(
    cl, varlist = c("spl")
  )
  
  for (i in 1:length(fusions)) {
    
    # find overlaps with fusion breakpoints:
    overlapping_reads <- parLapply(
      cl,
      spl,
      find_overlapping_reads,
      fusion = fusions[i],
      chromosome = "chr22"
    )
    
    overlapping_reads <- unlist(
      as(
        overlapping_reads[
          sapply(overlapping_reads, function(x) !is.null(x))
        ],
        "GRangesList"
      )
    )
    names(overlapping_reads) <- NULL
    
    if (!is.null(overlapping_reads)) {
      all_reads[[i]]$non_supporting$overlapping <- overlapping_reads
    } else {
      all_reads[[i]]$non_supporting$overlapping <- GRanges(NULL)
    }

  }
  
  stopCluster(cl)
  
  # find overlaps of split primary reads with chr22 breakpoint coords, 
  # add to supporting_$overlapping:
  spl <- split(filtered_reads$split_pairs, filtered_reads$split_pairs$qname)
  
  # initiate cluster:
  cl <- makeCluster(7)
  clusterExport(
    cl, varlist = c("spl")
  )
  
  for (i in 1:length(fusions)) {
    
    # find overlaps with fusion breakpoints:
    overlapping_reads <- parLapply(
      cl,
      spl,
      find_overlapping_reads,
      fusion = fusions[i],
      chromosome = "chr22"
    )
    
    overlapping_reads <- unlist(
      as(
        overlapping_reads[
          sapply(overlapping_reads, function(x) !is.null(x))
        ],
        "GRangesList"
      )
    )
    names(overlapping_reads) <- NULL
    
    if (!is.null(overlapping_reads)) {
      all_reads[[i]]$supporting$overlapping <- overlapping_reads
    } else {
      all_reads[[i]]$supporting$overlapping <- GRanges(NULL)
    }
    
  }
  
  # find overlaps of split primary reads with chr11 breakpoint coords, 
  # add to supporting_reads$overlapping:
  
  # find overlaps with fusion breakpoint:
  overlapping_reads <- parLapply(
    cl,
    spl,
    find_overlapping_reads,
    fusion = fusions[i],
    chromosome = "chr11"
  )
  
  stopCluster(cl)
  
  all_reads[[i]]$supporting$overlapping <- c(
    all_reads[[i]]$supporting$overlapping,
    unlist(
      as(
        overlapping_reads[
          sapply(overlapping_reads, function(x) !is.null(x))
        ],
        "GRangesList"
      )
    )
  )
  names(all_reads[[i]]$supporting$overlapping) <- NULL
  
  # filter reads:
  all_reads <- lapply(all_reads, function(x) {
    
    return(
      lapply(x, function(y) {
        
        if (length(x[[i]]$overlapping > 0)) {
        
          # deduplicate reads:
          spl <- split(y$overlapping, y$overlapping$qname)
          spl <- spl[!duplicated(spl)]
          
          y$overlapping <- unlist(
            as(
              spl[sapply(spl, function(z) !is.null(z))],
              "GRangesList"
            )
          )
          names(y$overlapping) <- NULL
          
          # filter overlapping reads without at least min_overlap on both sides of fusion:
          y$overlapping <- filter_overlaps(
            reads = y$overlapping, 
            min_overlap = min_overlap
          )
        }
        
        return(y)
        
      })
    )
    
  })

  
  ####################################################################################
  ### 3. Find breakpoint-spanning reads ###
  ####################################################################################
  
  # find breakpoint-spanning reads from non-split non-discordant reads, 
  # add to spanning$non_supporting:
  
  # fetch read gaps (read + gap coords):
  # split by qname:
  spl <- split(
    filtered_reads$non_split_concordant_pairs, 
    filtered_reads$non_split_concordant_pairs$qname
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
  
  # split by qname:
  spl <- split(non_split_concordant_gaps, non_split_concordant_gaps$qname)
  
  if (length(spl) > 0) {
    
    # initiate cluster:
    cl <- makeCluster(7)
    clusterExport(
      cl, varlist = c("spl")
    )
    
    for (i in 1:length(fusions)) {
      
      # find overlaps of non-supporting gaps with fusions:
      spanning_gaps <- parLapply(
        cl, 
        spl, 
        find_overlapping_reads, 
        fusions = fusions[i],
        chromosome = "chr22"
      )
      
      # merge ranges:
      spanning_gaps <- unlist(
        as(
          spanning_gaps[
            sapply(spanning_gaps, function(x) !is.null(x))
          ],
          "GRangesList"
        )
      )
      names(spanning_gaps) <- NULL
      
      stopCluster(cl)
      
      if (!is.null(spanning_gaps)) {
        
        # fetch corresponding read pairs:
        all_reads[[i]]$non_supporting$spanning <- filtered_reads$non_split_concordant_pairs[
          filtered_reads$non_split_concordant_pairs$qname %in% 
            spanning_gaps$qname
        ]
        m <- match(
          all_reads[[i]]$non_supporting$spanning$qname, spanning_gaps$qname
        )
        all_reads[[i]]$non_supporting$spanning$fusion_coord <- 
          spanning_gaps$fusion_coord[m]
        
      } else {
        all_reads[[i]]$non_supporting$spanning <- GRanges(NULL)
      }
      
    }
    
  } else {
    for (i in 1:length(fusions)) {
      all_reads[[i]]$non_supporting$spanning <- GRanges(NULL)
    }
  }
  
  # find breakpoint-spanning reads from discordant reads, 
  # add to spanning$supporting:
  spl <- split(
    filtered_reads$discordant_pairs, 
    filtered_reads$discordant_pairs$qname
  )
  
  for (i in 1:length(fusions)) {
    
    spanning_reads <- lapply(
      spl, 
      find_spanning_discordant, 
      fusions = fusions[i], 
      exp_window = disc_read_window
    )
    
    # merge to granges object:
    spanning_reads <- unlist(
      as(
        spanning_reads[
          sapply(spanning_reads, function(x) !is.null(x))
        ], "GRangesList"
      )
    )
    names(spanning_reads) <- NULL
    
    if (!is.null(overlapping_reads)) {
      all_reads[[i]]$supporting$spanning <- spanning_reads
    } else {
      all_reads[[i]]$supporting$spanning <- GRanges(NULL)
    }
    
  }
  
  
  ####################################################################################
  ### 4. Calculate proportions of non-supporting vs supporting read pairs for each
  # sample ###
  ####################################################################################
  
  # combine all supporting reads, and all non-supporting reads, for each fusion:
  combined_reads <- lapply(all_reads, function(x) {
    return(
      lapply(x, function(y) {
        return(c(y$spanning, y$overlapping))
      })
    )
  })
  
  # check all reads are in pairs:
  pair_check <- lapply(combined_reads, function(x) {
    sapply(x, function(y) {
      spl <- split(y, y$qname)
      return(
        all(
          sapply(spl, function(z) {
            length(z) == 2
          })
        )
      )
    })
  })
    
  print(
    paste0(
      "Are all supporting reads in pairs? ", 
      all(sapply(pair_check, function(x) x[names(x) == "supporting"]))
    )
  )
  
  print(
    paste0(
      "Are all non-supporting reads in pairs? ", 
      all(sapply(pair_check, function(x) x[names(x) == "non_supporting"]))
    )
  )
  
  save.image(paste0(Robject_dir, "temp_image.Rdata"))
  
  # save each group as sam file:
  for (i in 1:length(all_reads)) {
    for (j in 1:length(all_reads[[i]])) {
      for (k in 1:length(all_reads[[i]][[j]])) {
        
        if (length(all_reads[[i]][[j]][[k]]) > 0) {
          
          # define sam cols:
          sam <- data.frame(
            qname = all_reads[[i]][[j]][[k]]$qname,
            flag = all_reads[[i]][[j]][[k]]$flag,
            rname = seqnames(all_reads[[i]][[j]][[k]]),
            pos = start(all_reads[[i]][[j]][[k]]),
            mapq = all_reads[[i]][[j]][[k]]$mapq,
            cigar = all_reads[[i]][[j]][[k]]$cigar,
            rnext = all_reads[[i]][[j]][[k]]$rnext,
            pnext = all_reads[[i]][[j]][[k]]$pnext,
            tlen = all_reads[[i]][[j]][[k]]$tlen,
            seq = all_reads[[i]][[j]][[k]]$seq,
            qual = all_reads[[i]][[j]][[k]]$qual
          )
          
          # write sam to tab-separated file:
          write.table(
            sam,
            paste0(
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads_temp.sam"
            ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
          )
          
          # add header:
          system(
            paste0(
              "samtools view -H ", 
              bam_path, samplename, "/", samplename, ".consensus.concordant.pairs.bam",
              " > ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sam"
            )
          )
          
          # add rest of sam:
          system(
            paste0(
              "cat ", 
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads_temp.sam", 
              " >> ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sam"
            )
          )
          
          # convert to bam:
          system(
            paste0(
              "samtools view -bh ", 
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sam",
              " > ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.bam"
            )
          )
          
          # sort:
          system(
            paste0(
              "samtools sort -o ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sorted.bam ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.bam"
            )
          )
          
          # index:
          system(
            paste0(
              "samtools index ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sorted.bam"
            )
          )
          
          # clean:
          system(
            paste0(
              "rm ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads_temp.sam ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sam ",
              fusion_dir, "fusion_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.bam"
            )
          )
          
        }
        
      }
    }
  }

  saveRDS(all_reads, paste0(Robject_dir, "high_conf_fusion_reads.Rdata"))
  
} else {
  
  all_reads <- readRDS(paste0(Robject_dir, "high_conf_fusion_reads.Rdata"))
  
}

# keep breakpoints with greatest number of supporting reads:
read_counts <- sapply(combined_reads, function(x) {
  return(
    sapply(x, function(y) {
      length(y)
    })
  )
})

final_reads <- all_reads[[which.max(read_counts["supporting",])]]
final_combined <- combined_reads[[which.max(read_counts["supporting",])]]

final_venns <- lapply(final_reads, function(x) {
  return(lapply(x, create_venn, venn_cols))
})

# calculate supporting vs non-supporting proportion:
VAF <- round(
  length(final_combined$supporting)/length(final_combined$non_supporting)*100,
  1
)

# save VAF:
write.table(
  VAF,
  paste0(table_dir, "fusion", names(final_reads), "_VAF.txt"),
  quote = F,
  row.names = F,
  col.names = F
)
