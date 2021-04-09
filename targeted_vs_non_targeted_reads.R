project_name <- "ewing_ctDNA"

# define/create directories:
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
genome_dir <- paste0(project_dir, "genome/")
primer_dir <- paste0(
  project_dir, 
  "refs/QIAseq_DNA_panel.CDHS-34925Z-409.summary/"
)
result_dir <- paste0(project_dir, "results/")
gsnap_dir <- paste0(result_dir, "/gsnap/")

table_dir <- paste0(result_dir, "read_nos/")
system(paste0("mkdir -p ", table_dir))

samples <- list(
  "409_007_DB62M_TAGGCATG-CTCTCTAT_L001", 
  "409_012_DBV4V_CGTACTAG-CTCTCTAT_L001"
)


################################################################################
### 0. Load packages and functions ###
################################################################################

library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)


################################################################################
### 1. Load primer and sample bam co-ordinates ###
################################################################################

# load primer ranges:
captured_150 <- import(
  paste0(
  	primer_dir, 
  	"QIAseq_DNA_panel.CDHS-34925Z-409.primers-150bp.bed"
  ), 
  format = "bed"
)

captured_250 <- import(
  paste0(
  	primer_dir, 
  	"QIAseq_DNA_panel.CDHS-34925Z-409.primers-250bp.bed"
  ), 
  format = "bed"
)

# combine and reduce:
captured <- reduce(
  c(
  	captured_150,
  	captured_250
  )
)

# select bam fields to load:
what <- c("rname", "strand", "pos", "qwidth")
param <- ScanBamParam(what = what)

bams <- lapply(samples, function(x) {

  # load bams:
  bam <- scanBam(
  	BamFile(
  	  paste0(gsnap_dir, x, "/", x, ".all.sorted.bam"),
    ),
    param=param
  )

  # remove NAs:
  bam <- lapply(bam[[1]], function(y) y[!is.na(y)])

  # convert to GRanges:
  bam_gr <- GRanges(
  	seqnames = Rle(bam$rname),
  	ranges = IRanges(start = bam$pos, width = bam$qwidth),
  	strand = Rle(bam$strand)
  )

  # remove non-standard chromosomes:
  bam_gr <- bam_gr[grep("K|M|J|G", seqnames(bam_gr), invert = T)]

})

names(bams) <- samples


################################################################################
### 2. Find and remove reads overlapping with capture regions ###
################################################################################

for (i in 1:length(bams)) {

  olaps <- findOverlaps(bams[[i]], captured)

  if (i==1) {

  	non_capture_gr <- list(
  	  bams[[i]][-queryHits(olaps)]
  	)

  	read_nos <- data.frame(
  	  row.names = names(bams)[i],
  	  in_captures = length(unique(queryHits(olaps))),
  	  prop_in = round(length(unique(queryHits(olaps)))/length(bams[[i]]), 2),
  	  outside_captures = length(non_capture_gr[[i]]),
  	  prop_out = round(length(non_capture_gr[[i]])/length(bams[[i]]), 2)
  	)

  } else {

  	non_capture_gr[[i]] <- bams[[i]][-queryHits(olaps)]

  	read_nos <- rbind(
  	  read_nos,
  	  data.frame(
   	    row.names = names(bams)[i],
   	    in_captures = length(unique(queryHits(olaps))),
   	    prop_in = round(length(unique(queryHits(olaps)))/length(bams[[i]]), 2),
   	    outside_captures = length(non_capture_gr[[i]]),
   	    prop_out = round(length(non_capture_gr[[i]])/length(bams[[i]]), 2)
   	  )
  	)

  }

}


################################################################################
### 3. Find and count reads overlapping with fusion partner ###
################################################################################

if (!file.exists(paste0(genome_dir, "FLI1.gtf"))) {

  gencode <- import(paste0(genome_dir, "gencode.v19.annotation.gtf.gz"))

  FLI1_gr <- gencode[
    gencode$gene_name == "FLI1" & gencode$type == "gene"
  ]
  ETV1_gr <- gencode[
    gencode$gene_name == "ETV1" & gencode$type == "gene"
  ]

  export(FLI1_gr, paste0(genome_dir, "FLI1.gtf"))
  export(ETV1_gr, paste0(genome_dir, "ETV1.gtf"))

} else {

  FLI1_gr <- import(paste0(genome_dir, "FLI1.gtf"))
  ETV_gr <- import(paste0(genome_dir, "ETV1.gtf"))

}

read_nos$in_fusion <- c(
  sum(countOverlaps(bams[[1]], FLI1_gr)),
  sum(countOverlaps(bams[[2]], ETV1_gr))
)

# calculate proportion reads in fusion vs total:
read_nos$prop_in_fusion <- round(
  read_nos$in_fusion/
  (read_nos$in_captures + read_nos$outside_captures),
  5
)

# calculate proportion reads in fusion vs non-capture:
read_nos$prop_in_fusion_vs_non_capture <- round(
  read_nos$in_fusion/read_nos$outside_captures,
  5
)

# write to table:
write.table(
  read_nos,
  paste0(table_dir, "in_vs_out_captures_vs_fusion_read_nos.txt"),
  quote = F,
  row.names = T,
  col.names = T
)





