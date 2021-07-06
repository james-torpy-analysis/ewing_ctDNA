find_EWSR1_FLI1_fusions <- function(
  breakpoints, 
  GOI_list
) {

  # load function and define GOI regions:
  find_fusions <- dget(paste0(func_dir, "find_fusions.R"))

  all_GOI <- c(
    GOI_list$EWSR1,
    GOI_list$FLI1,
    GOI_list$ERG,
    GOI_list$ETV1,
    GOI_list$ETV4,
    GOI_list$FEV
  )

  # find breakpoint overlaps with EWSR1:
  EWSR1_fusions <- find_fusions(breakpoints, GOI_list$EWSR1)
  
  # count ranges:
  print("Number of EWSR1 fusions:")
  if (is.na(EWSR1_fusions)) {
    print(0)
  } else {
    print(length(EWSR1_fusions))
  }
  
  # determine whether joining ranges overlap FLI1, ETV1 or ERG:
  if (length(EWSR1_fusions) > 0 & !is.na(EWSR1_fusions)) {

    seqnams <- gsub(
      ":.*$", "", 
      gsub("^.*chr", "chr", EWSR1_fusions$join)
    )
    coord <- as.numeric(
      gsub(
        "[^0-9.-]", "", 
        gsub("^.*chr.*:", "", EWSR1_fusions$join)
      )
    )
    
    join_ranges <- GRanges(
      seqnames = EWSR1_fusions$join_chr,
      ranges = IRanges(
        start = EWSR1_fusions$join_coord, 
        end = EWSR1_fusions$join_coord
      ),
      strand = "*",
      join_chr = "chr22",
      join_coord = start(EWSR1_fusions)
    )

  } else {
    join_ranges <- NA
  }
  
  # identify EWSR1 fusions with GOI: 
  all_fusions <- list(
    true_positives = list(
      fusions = lapply(GOI_list, function(x) {
        return(find_fusions(join_ranges, x))
      })
    )
  )
    
  # count GOI fusions:
  all_fusions$true_positives$fusion_nos <- 
    lapply(all_fusions$true_positives$fusions, function(x) {
      if (all(is.na(x))) {
        return(NA)
      } else {
        return(length(x))
      }
    })

  # identify false EWSR1 fusions: 
  all_fusions$false_positives$fusions <- lapply(GOI_list, function(x) {
    return(find_fusions(join_ranges, x, invert = T))
  })
  
  # count GOI fusions:
  all_fusions$false_positives$fusion_nos <- 
    lapply(all_fusions$false_positives$fusions, function(x) {
      if (all(is.na(x))) {
        return(NA)
      } else {
        return(length(x))
      }
    })
   
  # clean up GOI_fusions:
  all_fusions <- lapply(all_fusions, function(x) {
    
    # change fusion NA values to empty GRanges:
    x$fusions[sapply(x$fusions, function(y) is.na(y)[1])] <- GRanges(NULL)
    
    # change fusion no NA values to 0:
    x$fusion_nos[sapply(x$fusion_nos, function(y) is.na(y)[1])] <- 0
      
    return(x)

  })
  
  return(all_fusions)

}
