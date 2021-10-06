longitudinal_patient_heatmap <- function(
  fusion_df, 
  hm_title,
  annot = "false positives",
  hm_cols
) {
  
  library(ComplexHeatmap)
  
  # change NAs in site column to 'unknown':
  fusion_df$Site[is.na(fusion_df$Site)] <- "unknown"
  
  # create VAF and detection dfs:
  VAF_df <- subset(fusion_df, select = c(Patient_id, Site, Treatment.dilution,
    Fusion_VAF ))
  
  detect_df <- VAF_df
  
  # make NA or not detected values blank:
  VAF_df$Fusion_VAF[
    is.na(VAF_df$Fusion_VAF) | VAF_df$Fusion_VAF == "not_detected"
  ] <- " "
  
  detect_df$Fusion_VAF[detect_df$Fusion_VAF != "not_detected"] <- "detected"

  VAF_df <- dcast(
    VAF_df, Patient_id + Site ~Treatment.dilution, value.var = "Fusion_VAF" )
  detect_df <- dcast(
    detect_df, Patient_id + Site ~Treatment.dilution, value.var = "Fusion_VAF" )
  
  # order dfs:
  VAF_df <- VAF_df[naturalorder(gsub("^.*_|M0|P", "", VAF_df$Patient_id)),]
  VAF_df$Site <- factor(VAF_df$Site, levels = c("primary", "met", "unknown"))
  VAF_df <- VAF_df[order(VAF_df$Site),]
  
  detect_df <- detect_df[naturalorder(gsub("^.*_|M0|P", "", detect_df$Patient_id)),]
  detect_df$Site <- factor(detect_df$Site, levels = c("primary", "met", "unknown"))
  detect_df <- detect_df[order(detect_df$Site),]
  
  # define site(row) and pathology(column) split vectors:
  met_split <- VAF_df$Site
  path_split <- factor(
    c("path", rep("non_path", length(unique(fusion_df$Treatment.dilution)))),
    levels = c("path", "non_path") )
  
  # make rownames are patient ids and add pathology column:
  rownames(VAF_df) <- VAF_df$Patient_id
  VAF_df$pathology <- " "
  VAF_df <- VAF_df[
    ,colnames(VAF_df) %in% 
      c("pathology", unique(as.character(fusion_df$Treatment.dilution))) ]
  VAF_df <- VAF_df %>%
    select(pathology, everything())
  
  # make rownames are patient ids and add pathology column:
  rownames(detect_df) <- detect_df$Patient_id
  detect_df <- subset(detect_df, select = -c(Patient_id, Site))
  path_df <- subset(fusion_df, select = c(Patient_id, Pathology_EWSR1_FLI1))
  path_df <- path_df[!duplicated(path_df$Patient_id),]
  path_df$Pathology_EWSR1_FLI1[
    path_df$Pathology_EWSR1_FLI1 == "yes" | grepl("chr", path_df$Pathology_EWSR1_FLI1)
  ] <- "pathology_detection"
  path_df$Pathology_EWSR1_FLI1[
    path_df$Pathology_EWSR1_FLI1 == "not_detected" ] <- "no_pathology_detection"
  path_df$Pathology_EWSR1_FLI1[is.na(path_df$Pathology_EWSR1_FLI1)] <- "unknown"
  rownames(path_df) <- path_df$Patient_id
  path_df <- subset(path_df, select = -Patient_id)
  colnames(path_df) <- "pathology"
  
  detect_df <- cbind(path_df, detect_df)
  
  # make NAs 'unknown':
  detect_df <- apply(detect_df, 2, function(x) {
    x[is.na(x)] <- "unknown"
    return(x)
  })
  
  # make NAs blank:
  VAF_df <- apply(VAF_df, 2, function(x) {
    x[is.na(x)] <- " "
    return(x)
  })
  
  # create heatmap:
  VAF_hm <- Heatmap(
    as.matrix(detect_df), 
    name = "Fusion detections", 
    row_split = met_split,
    column_split = path_split,
    col = hm_cols,
    border = "black",
    rect_gp = gpar(col = "black", lwd = 1),
    column_title = hm_title,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        VAF_df[i, j], x, y, 
        gp = gpar(fontsize = 10, fontface = "bold", col = "#991425") )
    }
  )
  
  # return heatmaps:
  return(list(VAF_hm, sread_hm))

}
