find_overlapping_reads <- function(
  read_pair,
  fusion,
  chromosome
) {
  
  # if finding overlaps with chr11 fusion coord, make this coord main one:
  if (chromosome == "chr11") {
    
    fusion <- GRanges(
      seqnames = fusion$join_chr,
      ranges = IRanges(start = fusion$join_coord, end = fusion$join_coord),
      strand = "*",
      join_chr = seqnames(fusion),
      join_coord = start(fusion)
    )
    
  }
  
  # find overlaps with fusion:
  olaps <- findOverlaps(read_pair, fusion)
  
  if (length(olaps) > 0) {
    
    # add overlapping breakpoint information:
    read_pair$fusion_chr <- NA
    read_pair$fusion_coord <- NA
    read_pair$fusion_mate_chr <- NA
    read_pair$fusion_mate_coord <- NA
    
    read_pair[queryHits(olaps)]$fusion_chr <- chromosome
    read_pair[queryHits(olaps)]$fusion_coord <- start(fusion)[
      subjectHits(olaps)
    ]
    
    read_pair[queryHits(olaps)]$fusion_mate_chr <- as.character(
      unique(fusion$join_chr)
    )
    read_pair[queryHits(olaps)]$fusion_mate_coord <- fusion$join_coord[
      subjectHits(olaps)
    ]
    
    # remove reads with end or start co-ordinates equal to breakpoint position,
    # as they do not span the breakpoints:
    if (
      all(
        start(read_pair[queryHits(olaps)]) == read_pair[queryHits(olaps)]$fusion_coord |
          end(read_pair[queryHits(olaps)]) == read_pair[queryHits(olaps)]$fusion_coord
      )
    ) {
      read_pair <- NULL
    }
    
  } else {
    read_pair <- NULL
  }
  
  return(read_pair)
  
}