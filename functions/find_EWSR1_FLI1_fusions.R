find_EWSR1_FLI1_fusions <- function(
  breakpoints, 
  GOI_list, 
  patient_df, 
  dilution_df
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
  EWSR1_fusions <- lapply(breakpoints, find_fusions, GOI_list$EWSR1)
  
  # remove NAs:
  EWSR1_fusions <- EWSR1_fusions[!is.na(EWSR1_fusions)]
  
  # count ranges:
  print("EWSR1 fusions per sample:")
  print(sapply(EWSR1_fusions, length))
  
  # determine whether joining ranges overlap FLI1, ETV1 or ERG:
  for (i in 1:length(EWSR1_fusions)) {

  
    if (length(EWSR1_fusions[[i]]) > 0) {
  
      seqnams <- gsub(
        ":.*$", "", 
        gsub("^.*chr", "chr", EWSR1_fusions[[i]]$join)
      )
      coord <- as.numeric(
        gsub(
          "[^0-9.-]", "", 
          gsub("^.*chr.*:", "", EWSR1_fusions[[i]]$join)
        )
      )
    
      if (!exists("join_ranges")) {
        join_ranges <- list(
          GRanges(
            seqnames = EWSR1_fusions[[i]]$join_chr,
            ranges = IRanges(
              start = EWSR1_fusions[[i]]$join_coord, 
              end = EWSR1_fusions[[i]]$join_coord
            ),
            strand = "*",
            join_chr = "chr22",
            join_coord = start(EWSR1_fusions[[i]])
          )
        )
      } else {
        join_ranges[[i]] <- GRanges(
          seqnames = EWSR1_fusions[[i]]$join_chr,
          ranges = IRanges(
            start = EWSR1_fusions[[i]]$join_coord, 
            end = EWSR1_fusions[[i]]$join_coord
          ),
          strand = "*",
          join_chr = "chr22",
          join_coord = start(EWSR1_fusions[[i]])
        )
      }
  
    } else {
  
      if (!exists("join_ranges")) {
        join_ranges <- list(NA)
      } else {
        join_ranges[[i]] <- NA
      }
  
    }  
  
  }
  names(join_ranges) <- names(EWSR1_fusions)
  
  # identify EWSR1 fusions with GOI: 
  GOI_fusions <- lapply(GOI_list, function(x) {
  
    for (i in 1:length(join_ranges)) {
    
      if (i==1) {
        fusions <- list(find_fusions(join_ranges[[i]], x))
      } else {
        fusions[[i]] <- find_fusions(join_ranges[[i]], x)
      }
    
    }
    names(fusions) <- names(join_ranges)
  
    return(fusions)
  
  })
  
  # count GOI fusions:
  fusion_nos <- lapply(GOI_fusions, function(x) {
  
    lengths <- sapply(x, function(y) {
      if (length(y) > 0) {
        if (!is.na(y)) {
          return(length(y))
        } else {
          return(0)
        }
      } else {
        return(0)
      }
    })
    names(lengths) <- gsub(
      "_.*$", "", 
      gsub("409_", "", names(x))
    )
  
    return(lengths)
  
  })

  # match fusion detections with sample numbers and add to patient_df:
  m <- match(patient_df$Sample, names(fusion_nos$FLI1))
  patient_df$Detected_FLI1_EWSR1_fusions = fusion_nos$FLI1[m]
  patient_df$Detected_FLI1_EWSR1_fusions[is.na(patient_df$Detected_FLI1_EWSR1_fusions)] <- 0

  # match fusion detections with sample numbers and add to patient_df:
  m <- match(dilution_df$Sample, names(fusion_nos$FLI1))
  dilution_df$Detected_FLI1_EWSR1_fusions = fusion_nos$FLI1[m]
  dilution_df$Detected_FLI1_EWSR1_fusions[is.na(dilution_df$Detected_FLI1_EWSR1_fusions)] <- 0

  # change EWSR1 fusions to yes/no if necessary:
  if (binary_calls) {
    patient_df$Detected_FLI1_EWSR1_fusions[as.numeric(patient_df$Detected_FLI1_EWSR1_fusions) == 0] <- "no"
    patient_df$Detected_FLI1_EWSR1_fusions[as.numeric(patient_df$Detected_FLI1_EWSR1_fusions) > 0] <- "yes"
    colnames(patient_df) <- gsub("fusions", "fusion", colnames(patient_df))
    dilution_df$Detected_FLI1_EWSR1_fusions[as.numeric(dilution_df$Detected_FLI1_EWSR1_fusions) == 0] <- "no"
    dilution_df$Detected_FLI1_EWSR1_fusions[as.numeric(dilution_df$Detected_FLI1_EWSR1_fusions) > 0] <- "yes"
    colnames(dilution_df) <- gsub("fusions", "fusion", colnames(dilution_df))
  }

  # identify false positive fusions: 
  for (i in 1:length(join_ranges)) {
  
    if (i==1) {
      false_fusions <- list(
        find_fusions(
          query_coord = join_ranges[[i]], 
          subject_coord = all_GOI, 
          invert = T
        )
      )
    } else {
      false_fusions[[i]] <- find_fusions(
        query_coord = join_ranges[[i]], 
        subject_coord = all_GOI, 
        invert = T
      )
    }
  
  }
  names(false_fusions) <- names(join_ranges)
  
  # count GOI fusions:
  false_fusion_nos <- sapply(false_fusions, function(y) {
    if (length(y) > 0) {
      if (!is.na(y)) {
        return(length(y))
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  })
  names(false_fusion_nos) <- gsub(
    "_.*$", "", 
    gsub("409_", "", names(false_fusions))
  )

  # match false fusion detections with sample numbers and add to patient_df:
  m <- match(patient_df$Sample, names(false_fusion_nos))
  patient_df$False_EWSR1_fusions = false_fusion_nos[m]
  patient_df$False_EWSR1_fusions[is.na(patient_df$False_EWSR1_fusions)] <- 0

  # match fusion detections with sample numbers and add to dilution_df:
  m <- match(dilution_df$Sample, names(false_fusion_nos))
  dilution_df$False_EWSR1_fusions = false_fusion_nos[m]
  dilution_df$False_EWSR1_fusions[is.na(dilution_df$False_EWSR1_fusions)] <- 0

  # clean up GOI_fusions:
  GOI_fusions <- GOI_fusions[!(names(GOI_fusions) %in% "EWSR1")]
  GOI_fusions <- lapply(GOI_fusions, function(x) {
    temp <- lapply(x, function(y) y[!is.na(y)])
    return(temp[lengths(temp) != 0])
  })
  return(
    list(
      fusion_nos = list(
        patient_df = patient_df,
        dilution_df = dilution_df
      ),
      GOI_fusions = GOI_fusions
    )
  )

}





