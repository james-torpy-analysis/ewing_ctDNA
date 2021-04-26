### find overlaps between gene SV breakpoints and other genes ###

find_fusions <- function(query_coord, subject_coord) {

  if (!is.na(query_coord)) {

    olaps <- findOverlaps(query_coord, subject_coord)

    if (length(olaps) > 0) {
      fusion <- query_coord[queryHits(olaps)]
    } else {
      fusion <- NA
    }

  } else {
    fusion <- NA
  }

  if (!is.na(fusion)) {
    # remove duplicates:
    fusion <- fusion[!duplicated(start(fusion))]
  
    # remove chr22 to chr22 fusions:
    fusion <- fusion[!( length(grep("chr22", seqnames(fusion))) > 0 & length(grep("chr22", fusion$join)) > 0 )]
  }

  return(fusion)
  
}