find_overlapping_reads <- function(
  read_pair,
  fusions,
  chromosome
) {
  
  # find overlaps with fusions:
  olaps <- findOverlaps(read_pair, fusions)
  
  if (length(olaps) > 0) {
    
    # add overlapping breakpoint information:
    read_pair$fusion_chr <- NA
    read_pair$fusion_coord <- NA
    read_pair$fusion_mate_chr <- NA
    read_pair$fusion_mate_coord <- NA
    
    read_pair[queryHits(olaps)]$fusion_chr <- chromosome
    read_pair[queryHits(olaps)]$fusion_coord <- start(fusions)[
      queryHits(olaps)
    ]
    
    read_pair[queryHits(olaps)]$fusion_mate_chr <- as.character(
      unique(fusions$join_chr)
    )
    read_pair[queryHits(olaps)]$fusion_mate_coord <- fusions$join_coord[
      queryHits(olaps)
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