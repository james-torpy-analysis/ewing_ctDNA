find_overlapping_reads <- function(
  read_pair,
  fusions,
  chromosome
) {
  
  # find overlaps with fusions:
  olaps <- findOverlaps(fusions, read_pair)
  
  if (length(olaps) > 0) {
    
    # add overlapping breakpoint information:
    read_pair$fusion_chr <- NA
    read_pair$fusion_coord <- NA
    read_pair[subjectHits(olaps)]$fusion_chr <- chromosome
    read_pair[subjectHits(olaps)]$fusion_coord <- start(fusions)[queryHits(olaps)]
    
    # remove reads with end or start co-ordinates equal to breakpoint position,
    # as they do not span the breakpoints:
    if (
      all(
        start(read_pair[subjectHits(olaps)]) == read_pair[subjectHits(olaps)]$fusion_coord |
          end(read_pair[subjectHits(olaps)]) == read_pair[subjectHits(olaps)]$fusion_coord
      )
    ) {
      read_pair <- NULL
    }
    
  } else {
    read_pair <- NULL
  }
  
  return(read_pair)
  
}