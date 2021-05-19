### find overlaps between gene SV breakpoints and other genes ###

find_fusions <- function(query_coord, subject_coord, invert = FALSE) {

  if (!is.na(query_coord[1])) {

    olaps <- findOverlaps(query_coord, subject_coord) 

    if (invert) {

      fusion <- query_coord[-queryHits(olaps)]
      if (length(fusion) < 1) {
        fusion <- NA
      }

    } else {

      if (length(olaps) > 0) {
        fusion <- query_coord[queryHits(olaps)]
      } else {
        fusion <- NA
      }

    }
    
  } else {
    fusion <- NA
  }

  if (!is.na(fusion)) {
    # remove duplicates:
    fusion <- fusion[!(duplicated(start(fusion)) & duplicated(fusion$join_coord))]
  }

  return(fusion)
  
}