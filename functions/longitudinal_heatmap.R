longitudinal_heatmap <- function(
  fusion_df, 
  hm_title, 
  type, 
  hm_cols,
  SNP_annot = FALSE
) {
  
  if (type == "patient") {

    # split df for each patient:
    split_df <- split(fusion_df, fusion_df$Patient)

    condition_vec <- c(
      "tumour", "naive", "NACT1", "NACT2", 
      "resection", "ACT1", "targeted", "relapse", 
      "ACT2"
    )

  } else if (type == "dilution") {

    # split df for each cell line:
    split_df <- split(fusion_df, fusion_df$Original_sample)

    condition_vec <- c(
      "100", "50", "10", "4", 
      "3", "2", "1", "0.4", 
      "0.2", "0.1", "0.04"
    )

  }
  
  # for each patient, make complete df of all treatments vs detections:
  hm_split <- lapply(split_df, function(x) {

    if (type == "patient") {

      empty_df <- data.frame(
        Treatment = factor(condition_vec, levels = condition_vec)
      )
      
      if (SNP_annot) {
        
        temp_sub <- subset(
          x, 
          select = c(
            Treatment, Detected_FLI1_EWSR1_fusion, 
            Known_EWSR1_FLI1_fusion, TP53_VAF, STAG2_VAF
          )
        )
        
      } else {
        
        temp_sub <- subset(
          x, 
          select = c(
            Treatment, Detected_FLI1_EWSR1_fusion, 
            Known_EWSR1_FLI1_fusion, False_EWSR1_fusions
          )
        )
        
      }
      
      merged_df <- merge(
        empty_df,
        temp_sub,
        by = "Treatment",
        all = T
      )
      merged_df <- merged_df %>%
        column_to_rownames("Treatment")

    } else if (type == "dilution") {

      empty_df <- data.frame(
        Dilution = factor(condition_vec, levels = condition_vec)
      )
      
      if (SNP_annot) {
        
        temp_sub <- subset(
          x, 
          select = c(
            Treatment, Detected_FLI1_EWSR1_fusion, 
            Known_EWSR1_FLI1_fusion, TP53_VAF, STAG2_VAF
          )
        )
        
      } else {
        
        temp_sub <- subset(
          x, 
          select = c(
            Treatment, Detected_FLI1_EWSR1_fusion, 
            Known_EWSR1_FLI1_fusion, False_EWSR1_fusions
          )
        )
        
      }

      merged_df <- merge(
        empty_df,
        temp_sub,
        by = "Dilution",
        all = T
      )

      # add 'b' to any duplicated dilutions:
      temp_dilution <- as.character(merged_df$Dilution)
      temp_dilution[duplicated(temp_dilution)] <- paste0(
        temp_dilution[duplicated(temp_dilution)], "b"
      )
      merged_df$Dilution <- temp_dilution
      merged_df <- merged_df %>%
        column_to_rownames("Dilution")

    }

    # order merged_df:
    m <- match(condition_vec, rownames(merged_df))
    merged_df <- merged_df[m,]
    
    if (SNP_annot) {
      
      # make zeros = NA:
      merged_df$TP53_VAF[merged_df$TP53_VAF == 0] <- NA
      merged_df$STAG2_VAF[merged_df$STAG2_VAF == 0] <- NA
      
      # merge SNP VAFs:
      merged_df$SNP_VAF <- NA
      
      for (j in 1:nrow(merged_df)) {
        
        if (!is.na(merged_df$TP53_VAF[j]) & !is.na(merged_df$STAG2_VAF[j])) {
          merged_df$SNP_VAF[j] <- paste0(
            "TP53_", merged_df$TP53_VAF[j], "_STAG2_", merged_df$STAG2_VAF[j]
          )
        } else if (!is.na(merged_df$TP53_VAF[j]) & is.na(merged_df$STAG2_VAF[j])) {
          merged_df$SNP_VAF[j] <- paste0("TP53_", merged_df$TP53_VAF[j])
        } else if (is.na(merged_df$TP53_VAF[j]) & !is.na(merged_df$STAG2_VAF[j])) {
          merged_df$SNP_VAF[j] <- paste0("STAG2_", merged_df$STAG2_VAF[j])
        }
        
      }
      
      # remove individual SNP_VAF columns:
      merged_df <- subset(merged_df, select = -c(TP53_VAF, STAG2_VAF))
      
    }

    # record FISH results for first sample:
    FISH = merged_df$Known_EWSR1_FLI1_fusion[
      !is.na(merged_df$Known_EWSR1_FLI1_fusion)
    ][1]
    # make FISH entries distinct from others:
    FISH <- paste0("FISH_", FISH)
    detection_df <- subset(
      merged_df, 
      select = Detected_FLI1_EWSR1_fusion
    )
    detection_df <- rbind(
      data.frame(row.names = "FISH", Detected_FLI1_EWSR1_fusion = FISH),
      detection_df
    )
    
    # create a SNP df in parallel:
    annot_df <- subset(
      merged_df, 
      select = -c(Known_EWSR1_FLI1_fusion, Detected_FLI1_EWSR1_fusion)
    )
    temp_bind <- as.data.frame(
      t(data.frame(rep(NA, ncol(annot_df))))
    )
    colnames(temp_bind) <- colnames(annot_df)
    rownames(temp_bind) <- "FISH"
    annot_df <- rbind(
      temp_bind,
      annot_df
    )
    
    return(list(detection_df = detection_df, annot_df = annot_df))
    
  })
  
  both_split <- list(
    detection_df = lapply(hm_split, function(x) x$detection_df),
    annot_df = lapply(hm_split, function(x) x$annot_df)
  )
  
  # bind together:
  hm_dfs <- lapply(both_split, function(x) {
    
    hm_df <- as.data.frame(t(do.call("cbind", x)))

    if (type == "patient") {

      # format data to order by met and FISH status:
      hm_df$Patient <- names(hm_split)
      met_status <- data.frame(
        Patient = fusion_df$Patient,
        Site = fusion_df$Site
      )
      met_status <- met_status[!duplicated(met_status$Patient),]
      hm_df <- merge(hm_df, met_status, by = "Patient")

      # split by primary/met:
      temp_split <- split(hm_df, hm_df$Site)

      if (!all(is.na(hm_df$FISH))) {
        temp_split <- lapply(temp_split, function(y) {
          # change all NAs to "unknown":
          y[is.na(y)] <- "unknown"
          # make FISH column a factor with required order:
          y$FISH <- factor(y$FISH, levels = c("FISH_yes", "FISH_no", "unknown"))
          # order each split by FISH column:
          y <- y[order(y$FISH),]
          return(y)
        })
      }
      names(temp_split) <- names(temp_split)

      # bind together in met status order:
      temp_split <- list(
        temp_split$primary,
        temp_split$met,
        temp_split$unknown
      )
      hm_df <- do.call("rbind", temp_split)
      # save local/met vector and move patient ids to rownames:
      met_order <- factor(hm_df$Site, levels = unique(hm_df$Site))
      hm_df <- subset(hm_df, select = -Site)
      rownames(hm_df) <- NULL
      hm_df <- hm_df %>%
        column_to_rownames("Patient")

      return(list(hm_df = hm_df, met_order = met_order))

    } else if (type == "dilution") {

      rownames(hm_df) <- paste0("Cell line sample ", names(hm_split))

      return(list(hm_df = hm_df))

    }
    
  })

  # adjust order of annot dataframe:
  hm_dfs$annot_df$hm_df <- hm_dfs$annot_df$hm_df[
    rownames(hm_dfs$detection_df$hm_df),
  ]
  
  # change all NA or 0 values in false positive df to spaces:
  final_annot <- apply(hm_dfs$annot_df$hm_df, 2, function(x) {
    x[is.na(x)] <- " "
    x[x == 0] <- " "
    return(x)
  })
  
  if (SNP_annot) {
    
    # remove underscores from SNP annotations:
    final_annot <- gsub(
      "_", "\n", 
      gsub("_STAG2", "\nSTAG2", final_annot)
    )
    
  }

  # create treatment_split vector:
  treatment_split <- c("FISH", rep("NOT_FISH", ncol(hm_dfs$detection_df$hm_df) - 1))
  
  # create heatmaps:
  if (type == "patient") {

    if (SNP_annot) {
      return(
        Heatmap(
          as.matrix(hm_dfs$detection_df$hm_df), 
          name = "Fusion detections", 
          row_split = hm_dfs$detection_df$met_order,
          column_split = treatment_split,
          col = hm_cols,
          border = "black",
          rect_gp = gpar(col = "black", lwd = 1),
          column_title = hm_title,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(
              final_annot[i, j], x, y, gp = gpar(
                fontsize = 5, fontface = "bold", col = "red"
              )
            )
          }
        )
      )
    } else {
      return(
        Heatmap(
          as.matrix(hm_dfs$detection_df$hm_df), 
          name = "Fusion detections", 
          row_split = hm_dfs$detection_df$met_order,
          column_split = treatment_split,
          col = hm_cols,
          border = "black",
          rect_gp = gpar(col = "black", lwd = 1),
          column_title = hm_title,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(
              final_annot[i, j], x, y, gp = gpar(fontsize = 10, col = "red")
            )
          }
        )
      )
    }
    

  } else if (type == "dilution") {

    if (SNP_annot) {
      return(
        Heatmap(
          as.matrix(hm_dfs$detection_df$hm_df), 
          name = "Fusion detections", 
          row_split = hm_dfs$detection_df$met_order,
          column_split = treatment_split,
          col = hm_cols,
          border = "black",
          rect_gp = gpar(col = "black", lwd = 1),
          column_title = hm_title,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(
              final_annot[i, j], x, y, gp = gpar(
                fontsize = 5, fontface = "bold", col = "red"
              )
            )
          }
        )
      )
    } else {
      return(
        Heatmap(
          as.matrix(hm_dfs$detection_df$hm_df), 
          name = "Fusion detections", 
          row_split = hm_dfs$detection_df$met_order,
          column_split = treatment_split,
          col = hm_cols,
          border = "black",
          rect_gp = gpar(col = "black", lwd = 1),
          column_title = hm_title,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(
              final_annot[i, j], x, y, gp = gpar(fontsize = 10, col = "red")
            )
          }
        )
      )
    }

  }
    
}