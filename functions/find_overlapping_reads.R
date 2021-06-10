find_overlapping_reads <- function(
  fusions,
  reads,
  chromosome
) {
  
  olaps <- findOverlaps(fusions, reads)
  
  overlapping <- reads[subjectHits(olaps)]
  
  # add overlapping breakpoint information:
  overlapping$fusion_coord <- start(fusions)[queryHits(olaps)]
  
  # remove reads with end or start co-ordinates equal to breakpoint position,
  # as they do not span the breakpoints:
  overlapping <- overlapping[
    start(overlapping) != overlapping$fusion_coord &
      end(overlapping) != overlapping$fusion_coord
  ]
  
  # remove duplicate read names, as only need one read of each pair:
  overlapping <- overlapping[
    !duplicated(overlapping$qname)
  ]
  
  # specify which chromosome breakpoint coordinates are on:
  colnames(mcols(overlapping))[
    colnames(mcols(overlapping)) == "fusion_coord"
  ] <- paste0(chromosome, "_fusion_coord")
  
  return(overlapping)
  
}