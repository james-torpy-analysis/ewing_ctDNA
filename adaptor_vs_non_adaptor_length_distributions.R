#!/bin/bash

sample_name <- "409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
fq_dir <- paste0(project_dir, "raw_files/temp/")
bam_path <- paste0(
  project_dir, "results/geneglobe/", sample_name, "/"
)
result_dir <- paste0(project_dir, "results/")

out_path <- paste0(result_dir, "QC/", sample_name, "/")
Robject_dir <- paste0(out_path, "Rdata/")
table_dir <- paste0(out_path, "tables/")
plot_dir <- paste0(out_path, "plots/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", plot_dir))


####################################################################################
### 0. Load packages ###
####################################################################################

library(Rsamtools)
library(ShortRead)
library(ggplot2)


####################################################################################
### 1. Load R1 fastq and bam ###
####################################################################################

# read in fastq:
fq_obj <- readFastq(paste0(fq_dir, sample_name, "_R1.fastq"))
fq <- data.frame(
  id = as.vector(id(fq_obj)),
  readseq = as.vector(sread(fq_obj))
)

# remove additional info from read names:
fq$id <- gsub(" .*$", "", fq$id)

# select fields to include in bam:
param <- ScanBamParam(
  what = c(
  	"qname", 
  	"flag",
  	"rname",
  	"seq",
  	"strand",
  	"pos",
  	"qwidth"
  )
)

# read in bam:
bam_obj <- scanBam(
  file = paste0(bam_path, sample_name, ".bam"), 
  index = paste0(bam_path, sample_name, ".bam"),
  param = param
)[[1]]
bam <- do.call("DataFrame", bam_obj)


####################################################################################
### 2. Plot overall read size distribution in bam ###
####################################################################################

# plot distribution of unfiltered mapped read lengths:
png(paste0(plot_dir, "mapped_read_length_hist.png"))
  plot(hist(width(bam$seq)))
dev.off()

# filter out reads with <40 bp length as recommended by Qiagen (primer dimers):
filt_bam <- bam[bam$qwidth >= 40,]

# plot distribution of filtered mapped read lengths:
png(paste0(plot_dir, "filtered_mapped_read_length_hist.png"))
  plot(hist(width(filt_bam$seq)))
dev.off()

filt_fq <- fq[fq$id %in% filt_bam$qname,]


####################################################################################
### 3. Identify adaptor vs non-adaptor reads from fastq ###
####################################################################################

adaptor_reads <- filt_fq[
  grep("^CCACCTGAAGTCCAAAAAGGGT|^GGATATATTTAAAGAAATGAATGACTTTTAACGAGA", filt_fq$readseq),
]

non_adaptor_reads <- filt_fq[
  grep(
    "^CCACCTGAAGTCCAAAAAGGGT|^GGATATATTTAAAGAAATGAATGACTTTTAACGAGA", 
    filt_fq$readseq,
    invert = TRUE
  ),
]

barplot_df <- data.frame(
  type = c("Adaptor", "Non-adaptor"),
  value = c(nrow(adaptor_reads), nrow(non_adaptor_reads)),
  x = c(1, 1)
)

p <- ggplot(barplot_df, aes(fill=type, y=value, x=x))
p <- p + geom_bar(position="dodge", stat="identity")
p <- p + ylab("log10 read no.")

png(paste0(plot_dir, "proportion_adaptor_vs_non_adaptor_reads.png"))
  p
dev.off()


####################################################################################
### 4. Divide bam into non-adaptor trimmed vs adaptor trimmed
# and plot distribution of lengths ###
####################################################################################



