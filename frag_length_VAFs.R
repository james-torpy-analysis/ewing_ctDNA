
projectname <- "ewing_ctDNA"
samplename <- "409_014_DBV4V_TCCTGAGC-CTCTCTAT_L001"

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
result_dir <- paste0(project_dir, "results/")
bam_dir <- paste0(result_dir, "BWA_and_picard/bams/")
out_path <- paste0(result_dir, "VAF_calculation/", samplename, "/")

Robject_dir <- paste0(out_path, "Rdata/")
plot_dir <- paste0(out_path, "plots/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

func_dir <- paste0(project_dir, "scripts/functions/")

library(GenomicAlignments)
library(Rsamtools)
# library(scales)
# 
# source(paste0(func_dir, "vaf_functions.R"))

# load filtered reads:
filtered_reads <- readRDS(paste0(Robject_dir, "filtered_reads.Rdata"))

# split by read name:
split_gr <- split(filtered_reads, filtered_reads$qname)

# remove entries with number of reads != 2:
length(split_gr)
split_gr <- split_gr[sapply(split_gr, function(x) length(x) == 2)]
length(split_gr)

# remove entries with reads mapping to any chr other than chr11 or chr22:
split_gr <- split_gr[sapply(split_gr, function(x) {
  return(all(seqnames(x) %in% c("chr11", "chr22")))
})]
length(split_gr)

# remove entries with either 2x + or 2x - reads:
split_gr <- split_gr[sapply(split_gr, function(x) {
  return(all(!duplicated(x$strand)))
})]
length(split_gr)

# keep only read pairs with both mapping to chr22:
chr22_split_gr <- split_gr[sapply(split_gr, function(x) {
  return(all(seqnames(x) == "chr22"))
})]
length(chr22_split_gr)

# calculate fragment lengths:
frag_lengths <- sapply(chr22_split_gr, function(x) length(start(x)[1]:end(x)[2]))
summary(frag_lengths)

# remove reads with length > 200:
frag_lengths <- frag_lengths[frag_lengths <= 200]
summary(frag_lengths)

plot(density(as.numeric(frag_lengths)))

png(paste0(plot_dir, "fragment_length_density.png"))
  plot(density(as.numeric(frag_lengths)))
dev.off()

# split into reads below 159 and those above:
normal <- length(frag_lengths[frag_lengths >= 159])
cancer <- length(frag_lengths[frag_lengths < 159])

