
sample_name <- "409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
raw_dir <- paste0(project_dir, "raw_files/")
fq_dir <- paste0(result_dir, "picard/fastqs/")
int_dir <- paste0(result_dir, "picard/int_bams/")
bam_dir <- paste0(result_dir, "picard/bams/")
result_dir <- paste0(project_dir, "results/")
stats_dir <- paste0(result_dir, "stats/", sample_name, "/")
temp_dir <- paste0(result_dir, "picard/temp/")

system(paste0("mkdir -p ", stats_dir))
system(paste0("mkdir -p ", temp_dir))


##################################################################################################################################
### 0. Load libraries and functions ###
##################################################################################################################################

library(ShortRead)
library(rtracklayer)
library(GenomicRanges)

mround <- function(x,base){
  base*round(x/base)
}


##################################################################################################################################
### 1. Load data ###
##################################################################################################################################

# load fastqs:
fqs <- list(
  orig_fq = readFastq(paste0(raw_dir, sample_name, "_R1.fastq.gz")),
  filtered_orig_fq = readFastq(paste0(fq_dir, sample_name, ".withUMI.fastq")),
  consensus_fq = readFastq(paste0(fq_dir, sample_name, ".consensus.fastq"))
)


mapped_bams <- list(
  initial_mapped_bam = import(paste0(bam_dir, sample_name, ".uncollapsed.bam")),
  initial_mapped_discordant_bam = import(paste0(bam_dir, sample_name, ".uncollapsed.discordant.bam")),
  initial_mapped_split_bam = import(paste0(bam_dir, sample_name, ".uncollapsed.split.bam")),
  consensus_remapped_bam = import(paste0(bam_dir, sample_name, ".consensus.unfiltered.bam")),
  consensus_remapped_discordant_bam = import(paste0(bam_dir, sample_name, ".consensus.unfiltered.discordant.bam")),
  consensus_remapped_split_bam = import(paste0(bam_dir, sample_name, ".consensus.unfiltered.split.bam")),
  filtered_consensus_bam = import(paste0(bam_dir, sample_name, ".consensus.filtered.bam")),
  filtered_consensus_discordant_bam = import(paste0(bam_dir, sample_name, ".consensus.filtered.discordant.bam")),
  filtered_consensus_split_bam = import(paste0(bam_dir, sample_name, ".consensus.filtered.split.bam"))
)


##################################################################################################################################
### 2. Determine sequence lengths and total read no from original fq1 file and collapsed fq file ###
##################################################################################################################################

fq_lengths <- lapply(fqs, function(x) {
  return(width(sread(x)))
})


##################################################################################################################################
### 3. Determine sequence lengths and total raw and read nos from bams ###
##################################################################################################################################

# determine unmapped bam names:
unmapped_names <- c(".unmapped.bam", ".unmapped.filtered.bam", ".consensus.unmapped.bam")
for (i in 1:length(unmapped_names)) {

  print(i)

  # remove header:
  system(
    paste0(
    	"samtools view ", int_dir, sample_name, unmapped_names[i], " > ", temp_dir, sample_name, "_no_header", 
      unmapped_names[i]
    )
  )
  # load bam:
  df <- read.table(
    paste0(temp_dir, sample_name, "_no_header", unmapped_names[i]),
    sep = "\t",
    fill = TRUE,
    header = FALSE
  )
  	
  # create length list:
  if (i==1) {
    unmapped_bam_lengths <- list(nchar(df$V10))
  } else {
    unmapped_bam_lengths[[i]] <- nchar(df$V10)
  }

}
names(unmapped_bam_lengths) <- gsub(
  "^_", "", 
  gsub("\\.", "_", unmapped_names)
)

mapped_bam_lengths <- lapply(mapped_bams, function(x) {
  return(qwidth(x))
})


##################################################################################################################################
### 4. Combine and plot lengths ###
##################################################################################################################################

# chronologically order outputs:
all_lengths <- list(
  orig_fq = fq_lengths$orig_fq,
  unmapped_bam = unmapped_bam_lengths$unmapped_bam,
  unmapped_filtered_bam = unmapped_bam_lengths$unmapped_filtered_bam,
  filtered_orig_fq = fq_lengths$filtered_orig_fq,
  initial_mapped_bam = mapped_bam_lengths$initial_mapped_bam,
  initial_mapped_discordant_bam = mapped_bam_lengths$initial_mapped_discordant_bam,
  initial_mapped_split_bam = mapped_bam_lengths$initial_mapped_split_bam,
  consensus_unmapped_bam = unmapped_bam_lengths$consensus_unmapped_bam,
  consensus_fq = fq_lengths$consensus_fq,
  consensus_remapped_bam = mapped_bam_lengths$consensus_remapped_bam,
  consensus_remapped_discordant_bam = mapped_bam_lengths$consensus_remapped_discordant_bam,
  consensus_remapped_split_bam = mapped_bam_lengths$consensus_remapped_split_bam,
  filtered_consensus_bam = mapped_bam_lengths$filtered_consensus_bam,
  filtered_consensus_discordant_bam = mapped_bam_lengths$filtered_consensus_discordant_bam,
  filtered_consensus_split_bam = mapped_bam_lengths$filtered_consensus_split_bam
)

for (i in 1:length(all_lengths)) {

  png(paste0(stats_dir, i, ".", names(all_lengths)[i], "_read_length_distribution.png"))
    plot(density(all_lengths[[i]]))
  dev.off()

}








