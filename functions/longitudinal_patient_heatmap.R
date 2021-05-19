longitudinal_patient_heatmap <- function(fusion_df, hm_title) {
  
  # split df for each patient:
  split_df <- split(fusion_df, fusion_df$Patient)
  treatment_vec <- c(
    "tumour", "resection", "naive", "NACT1", 
    "NACT2", "targeted", "ACT1", "relapse", 
    "ACT2"
  )
  # for each patient, make complete df of all treatments vs detections:
  hm_split <- lapply(split_df, function(x) {
    
    empty_df <- data.frame(
      Treatment = factor(treatment_vec, levels = treatment_vec)
    )
    merged_df <- merge(
      empty_df,
      subset(
        x, 
        select = c(
          Treatment, Detected_FLI1_EWSR1_fusion, 
          Known_EWSR1_FLI1_fusion, False_EWSR1_fusions
        )
      ),
      by = "Treatment",
      all = T
    )
    merged_df <- merged_df %>%
      column_to_rownames("Treatment")
    
    # record fish results for first sample:
    fish = merged_df$Known_EWSR1_FLI1_fusion[
      !is.na(merged_df$Known_EWSR1_FLI1_fusion)
    ][1]
    detection_df <- subset(
      merged_df, 
      select = -c(Known_EWSR1_FLI1_fusion, False_EWSR1_fusions)
    )
    detection_df <- rbind(
      data.frame(row.names = "fish", Detected_FLI1_EWSR1_fusion = fish),
      detection_df
    )
    
    # create a false positive df in parallel:
    fp_df <- subset(
      merged_df, 
      select = -c(Known_EWSR1_FLI1_fusion, Detected_FLI1_EWSR1_fusion)
    )
    fp_df <- rbind(
      data.frame(row.names = "fish", False_EWSR1_fusions = NA),
      fp_df
    )
    
    return(list(detection_df = detection_df, fp_df = fp_df))
    
  })
  both_split <- list(
    detection_df = lapply(hm_split, function(x) x$detection_df),
    fp_df = lapply(hm_split, function(x) x$fp_df)
  )
  
  # bind together:
  hm_dfs <- lapply(both_split, function(x) {
    
    hm_df <- as.data.frame(t(do.call("cbind", x)))
    hm_df$Patient <- names(hm_split)
    
    # group by primary/met:
    met_status <- data.frame(
      Patient = fusion_df$Patient,
      Site = fusion_df$Site
    )
    met_status <- met_status[!duplicated(met_status$Patient),]
    hm_df <- merge(hm_df, met_status, by = "Patient")
    temp_split <- split(hm_df, hm_df$Site)
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
    
  })
  
  # change all NA or 0 values in false positive df to spaces:
  final_fp <- apply(hm_dfs$fp_df$hm_df, 2, function(x) {
    x[is.na(x)] <- " "
    x[x == 0] <- " "
    return(x)
  })
  
  # create heatmaps:
  return(
    Heatmap(
      as.matrix(hm_dfs$detection_df$hm_df), 
      name = "fusion_detections", 
      row_split = hm_dfs$detection_df$met_order,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = hm_title,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          final_fp[i, j], x, y, gp = gpar(fontsize = 10, col = "red")
        )
      }
    )
  )
  
}