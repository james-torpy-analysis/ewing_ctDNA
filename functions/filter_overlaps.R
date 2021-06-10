filter_overlaps <- function(reads, min_overlap) {
  
  # split by chromosome:
  split_reads <- split(reads, seqnames(reads))
  split_reads <- split_reads[c("chr11", "chr22")]
  
  # calculate lengths from start of read to fusion, and from fusion to end:
  split_reads <- lapply(split_reads, function(y) {
    if (length(y) > 0) {
      if (unique(seqnames(y)) == "chr11") {
        y$start_to_fusion <- y$chr11_fusion_coord - start(y)
        y$fusion_to_end <- end(y) - y$chr11_fusion_coord
      } else if (unique(seqnames(y)) == "chr22") {
        y$start_to_fusion <- y$chr22_fusion_coord - start(y)
        y$fusion_to_end <- end(y) - y$chr22_fusion_coord
      }
    }
    return(y)
  })
  
  # concatentate reads back together:
  prefilt_reads <- c(split_reads[[1]], split_reads[[2]])
  
  # filter out reads without overlaps of at least min_overlap on both sides of
  # fusion:
  return(
    prefilt_reads[
      prefilt_reads$start_to_fusion >= min_overlap & 
        prefilt_reads$fusion_to_end >= min_overlap
    ]
  )
  
}