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