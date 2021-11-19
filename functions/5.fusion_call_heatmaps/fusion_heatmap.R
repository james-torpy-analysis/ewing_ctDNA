fusion_heatmap <- function(
  fusion_df,
  type,
  hm_cols
) {
  
  library(naturalsort)
  library(reshape2)
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

  # cast to wide format:
  VAF_df <- dcast(
    VAF_df, Patient_id + Site ~Treatment.dilution, value.var = "Fusion_VAF" )
  detect_df <- dcast(
    detect_df, Patient_id + Site ~Treatment.dilution, value.var = "Fusion_VAF" )
  
  # make rownames patient ids and add pathology column:
  rownames(VAF_df) <- VAF_df$Patient_id
  VAF_df$pathology <- " "
  VAF_df <- VAF_df[
    ,colnames(VAF_df) %in% 
      c("pathology", unique(as.character(fusion_df$Treatment.dilution))) ]
  VAF_df <- VAF_df %>%
    select(pathology, everything())
  
  # merge pathology results with dfs:
  rownames(detect_df) <- detect_df$Patient_id
  detect_df <- subset(detect_df, select = -Patient_id)
  path_df <- subset(fusion_df, select = c(Patient_id, Pathology_EWSR1_FLI1))
  path_df <- path_df[!duplicated(path_df$Patient_id),]
  path_df$Pathology_EWSR1_FLI1[
    path_df$Pathology_EWSR1_FLI1 == "yes" | grepl("chr", path_df$Pathology_EWSR1_FLI1)
  ] <- "pathology_detection"
  path_df$Pathology_EWSR1_FLI1[
    path_df$Pathology_EWSR1_FLI1 == "not_detected" ] <- "no_pathology_detection"
  path_df$Pathology_EWSR1_FLI1[is.na(path_df$Pathology_EWSR1_FLI1)] <- "no_sample"
  rownames(path_df) <- path_df$Patient_id
  path_df <- path_df[match(rownames(detect_df), path_df$Patient_id),]
  path_df <- subset(path_df, select = -Patient_id)
  colnames(path_df) <- "pathology"
  
  detect_df <- cbind(path_df, detect_df)
  
  # order dfs by name, pathology status, and site:
  if (type == "patient") {
    detect_df <- as.data.frame(
      detect_df[naturalorder(gsub("^.*_|M0|P", "", rownames(detect_df))),] )
    detect_df$pathology <- factor(detect_df$pathology, 
                                  levels = c("pathology_detection", "no_pathology_detection", "no_sample") ) 
    detect_df <- detect_df[order(detect_df$pathology),]
    detect_df$Site <- factor(
      detect_df$Site, levels = c("primary", "met", "unknown") )
    detect_df <- detect_df[order(detect_df$Site),]
  } else {
    detect_df <- detect_df[order(rownames(detect_df)),]
  }
  
  # define site(row) and pathology(column) split vectors:
  if (type == "patient") {
    met_split <- detect_df$Site
  }
  path_split <- factor(
    c("path", rep("non_path", length(unique(fusion_df$Treatment.dilution)))),
    levels = c("path", "non_path") )
  detect_df <- subset(detect_df, select = -Site)
  
  VAF_df <- VAF_df[match(rownames(detect_df), rownames(VAF_df)),]
  
  # make NAs 'no_sample':
  detect_df <- apply(detect_df, 2, function(x) {
    x[is.na(x)] <- "no_sample"
    return(x)
  })
  
  # remove underscores:
  rownames(detect_df) <- gsub("_", " ", rownames(detect_df))
  detect_df <- gsub("_", " ", detect_df)
  names(hm_cols) <- gsub("_", " ", names(hm_cols))
  
  # create annot_df:
  annot_df <- apply(VAF_df, 2, function(x) {
    x[!is.na(x) & x != " "] <- paste0("EWSR1/FLI1: ", x[!is.na(x) & x != " "])
    return(x)
  })
  
  # make NAs blank:
  annot_df <- apply(annot_df, 2, function(x) {
    x[is.na(x)] <- " "
    return(x)
  })
  
  if (type == "dilution") {
    # order by dilution level:
    detect_df <- detect_df[, 
      c("pathology", unique(fusion_df$Treatment.dilution)) ]
    annot_df <- annot_df[,
      c("pathology", unique(fusion_df$Treatment.dilution)) ]
  }
  
  # create VAF heatmap:
  if (type == "patient") {
    VAF_hm <- Heatmap(
      as.matrix(detect_df), 
      name = "Fusion detections", 
      row_split = met_split,
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient EWSR1/FLI1 fusion detections",
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          annot_df[i, j], x, y, 
          gp = gpar(fontsize = 6, fontface = "bold", col = "#991425") )
      }, )
  } else {
    VAF_hm <- Heatmap(
      as.matrix(detect_df), 
      name = "Fusion detections",
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient EWSR1/FLI1 fusion detections",
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          annot_df[i, j], x, y, 
          gp = gpar(fontsize = 6, fontface = "bold", col = "#991425") )
      }, )
  }
  
  # create supporting read annotation:
  supp_df <- subset(fusion_df, select = c(Patient_id, Site, Treatment.dilution,
    Forward_supporting ))
  
  # make NA or not detected values blank:
  supp_df$Forward_supporting[
    is.na(supp_df$Forward_supporting) | supp_df$Forward_supporting == "not_detected"
  ] <- " "
  
  # cast to wide format:
  supp_df <- dcast(
    supp_df, Patient_id + Site ~Treatment.dilution, value.var = "Forward_supporting" )
  
  # make rownames patient ids and add pathology column:
  rownames(supp_df) <- supp_df$Patient_id
  supp_df$pathology <- " "
  supp_df <- supp_df[
    ,colnames(supp_df) %in% 
      c("pathology", unique(as.character(fusion_df$Treatment.dilution))) ]
  supp_df <- supp_df %>%
    select(pathology, everything())
  
  # remove underscores:
  rownames(supp_df) <- gsub("_", " ", rownames(supp_df))
  
  # order df:
  supp_df <- supp_df[match(rownames(detect_df), rownames(supp_df)), 
    match(colnames(detect_df), colnames(supp_df)) ]
  
  # make NAs blank:
  supp_df <- apply(supp_df, 2, function(x) {
    x[is.na(x)] <- " "
    return(x)
  })
  
  # create VAF heatmap:
  if (type == "patient") {
    sread_hm <- Heatmap(
      as.matrix(detect_df), 
      name = "Fusion-supporting reads", 
      row_split = met_split,
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient EWSR1/FLI1 fusion-supporting reads",
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          supp_df[i, j], x, y, 
          gp = gpar(fontsize = 8, fontface = "bold", col = "#430F82") )
      },
    )
  } else {
    sread_hm <- Heatmap(
      as.matrix(detect_df), 
      name = "Fusion-supporting reads", 
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient EWSR1/FLI1 fusion-supporting reads",
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          supp_df[i, j], x, y, 
          gp = gpar(fontsize = 8, fontface = "bold", col = "#430F82") )
      },
    )
  }

  # return heatmaps:
  return(list(VAF=VAF_hm, sread=sread_hm))

}
