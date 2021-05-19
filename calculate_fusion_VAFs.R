
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
func_dir <- paste0(project_dir, "scripts/functions/")

fusion_dir <- paste0(project_dir, "results/fusions/")
bam_path <- paste0(project_dir, "results/BWA_and_picard/bams/")

####################################################################################
### 0. Load packages and functions ###
####################################################################################

library(rtracklayer)
library(GenomicRanges)


####################################################################################
### 1. Load and format data ###
####################################################################################

# load in fusions:
fusions <- readRDS(paste0(fusion_dir, "EWSR1_GOI_fusions.Rdata"))
fusions <- fusions$FLI1

# load in bam, split bam and discordant bam as GRanges:
samplenames <- names(fusions)
bams <- lapply(samplenames, function(x) {

  # define file types to import:
  bamtypes <- list(
  	all = "consensus.bam", 
  	split = "consensus.split.bam", 
  	discordant = "consensus.discordant.bam"
  )

  # import bams:
  bam_qnames <- lapply(bamtypes, function(y) {
  	return( scanBam(paste0(bam_path, x, "/", x, ".", y)) )
  })

  # add fields to granges:
  

})
names(bams) <- samplenames

# create split_primary_alignment GRanges:


# create non_split_or_discordant GRanges:

# create empty list to be filled:
#[[1]] supporting
#[[1]] [[1]] spanning_reads
#[[1]] [[2]] overlapping_reads
#[[2]] non_supporting
#[[2]] [[1]] spanning_reads
#[[2]] [[2]] overlapping_reads


####################################################################################
### 2. find breakpoint-overlapping reads: ###
####################################################################################

# create read list:
#[[1]] non_supporting
#[[2]] supporting_split

# find overlaps of reads with breakpoints, add to each 'overlapping_reads' element:

# for overlapping reads, calculate overlaps on each side:

# keep reads with breakpoint - start_coordinate >= 19 bp

# keep reads with end_coordinate - breakpoint >= 19 bp


####################################################################################
### 3. Find breakpoint-spanning reads ###
####################################################################################

# create non_supporting_fragments (from non_split_or_discordant):

# find overlaps of non_supporting_fragments and add to 'spanning_reads' element:

# find discordant reads which overlap breakpoint, and add to 'spanning reads' 
# element:
# keep reads within 300 bp of chr22 breakpoint:

# keep reads within 300 bp of chr11 breakpoint:


####################################################################################
### 4. Calculate proportions of non-supporting vs supporting read pairs for each
# sample ###
####################################################################################

# keep breakpoint with greatest number of supporting reads:

# add both supporting and non-supporting filtered overlapping and spanning reads


# calculate supporting vs non-supporting proportion

