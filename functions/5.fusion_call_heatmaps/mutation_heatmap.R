mutation_heatmap <- function(
  fusion_df,
  type,
  hm_cols
) {
  
  library(naturalsort)
  library(reshape2)
  library(ComplexHeatmap)
  
  round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5 + sqrt(.Machine$double.eps)
    z = trunc(z)
    z = z/10^n
    z*posneg
  }
  
  # change NAs in site column to 'unknown':
  fusion_df$Site[is.na(fusion_df$Site)] <- "unknown"
  
  # create VAF and detection dfs:
  VAF_df <- subset(fusion_df, select = c(Patient_id, Site, Treatment.dilution))
  VAF_df$mutation_VAF <- apply(fusion_df, 1, function(x) {
    x <- x[grepl("smCounter2|Fusion", names(x)) & grepl("VAF", names(x))]
    x <- x[grep("not_detected", x, invert = TRUE)]
    
    if (length(x) > 0) {
      # round VAFs:
      VAF <- round2(as.numeric(x), 1)
      names(VAF) <- names(x)
      names(VAF) <- gsub("Fusion", "EWSR1/FLI1", names(VAF))
      for (i in seq_along(VAF)) {
        if (i==1) {
          annot <- paste0(gsub(
            "smCounter2_|_VAF", "", names(VAF)[i]), ": ", VAF[i], "%" )
        } else {
          annot <- paste0(annot, 
            "\n", gsub("smCounter2_|_VAF", "", names(VAF)[i]), ": ", VAF[i], "%" )
        }
      }
      return(annot)
    } else {
      return("")
    }
    
  })
  
  detect_df <- VAF_df
  
  # annotate detected/not_detected values:
  detect_df$mutation_VAF[detect_df$mutation_VAF == ""] <- "not_detected"
  detect_df$mutation_VAF[detect_df$mutation_VAF != "not_detected"] <- "detected"

  # cast to wide format:
  VAF_df <- dcast(
    VAF_df, Patient_id + Site ~Treatment.dilution, value.var = "mutation_VAF" )
  detect_df <- dcast(
    detect_df, Patient_id + Site ~Treatment.dilution, value.var = "mutation_VAF" )
  
  # make rownames patient ids and add pathology column:
  rownames(VAF_df) <- VAF_df$Patient_id
  VAF_df$pathology <- ""
  VAF_df <- VAF_df[
    ,colnames(VAF_df) %in% 
      c("pathology", unique(as.character(fusion_df$Treatment.dilution))) ]
  VAF_df <- VAF_df %>%
    select(pathology, everything())
  
  # merge pathology results with dfs:
  rownames(detect_df) <- detect_df$Patient_id
  detect_df <- subset(detect_df, select = -Patient_id)
  path_df <- subset(fusion_df, select = c(
    Patient_id, Sanger_TP53_point_mut, 
    Sanger_STAG2_point_mut, Pathology_EWSR1_FLI1 ))
  path_df <- path_df[!duplicated(path_df$Patient_id),]
  
  # if were are any mutations picked up, indicate "pathology_detection":
  path_df$pathology <- "no_pathology_detection"
  path_df$pathology[apply(
    subset(path_df, select = -c(Patient_id, pathology)), 1, function(x) {
      any(x != "not_detected")
    })] <- "pathology_detection"
  path_df$pathology[apply(
    subset(path_df, select = -c(Patient_id, pathology)), 1, function(x) {
      all(is.na(x))
    })] <- "no_sample"
  
  # order as in detect_df:
  path_df <- path_df[match(rownames(detect_df), path_df$Patient_id),]
  rownames(path_df) <- path_df$Patient_id
  path_df <- subset(path_df, select = pathology)
  
  # bind to detect_df:
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
  
  # make NAs blank:
  VAF_df <- apply(VAF_df, 2, function(x) {
    x[is.na(x)] <- ""
    return(x)
  })
  
  if (type == "dilution") {
    # order by dilution level:
    detect_df <- detect_df[, 
      c("pathology", unique(fusion_df$Treatment.dilution)) ]
    VAF_df <- VAF_df[,
      c("pathology", unique(fusion_df$Treatment.dilution)) ]
    # add percentages to dilutions:
    colnames(detect_df)[colnames(detect_df) != "pathology"] <- paste0(
      colnames(detect_df)[colnames(detect_df) != "pathology"], "%" )
    colnames(VAF_df)[colnames(VAF_df) != "pathology"] <- paste0(
      colnames(VAF_df)[colnames(VAF_df) != "pathology"], "%" )
  }
  
  # create VAF heatmap:
  if (type == "patient") {
    VAF_hm <- Heatmap(
      as.matrix(detect_df), 
      name = "Mutation detections", 
      row_split = met_split,
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient mutation detections",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          VAF_df[i, j], x, y, 
          gp = gpar(fontsize = 7.5, fontface = "bold", col = "#991425") )
      }, )
  } else {
    VAF_hm <- Heatmap(
      as.matrix(detect_df), 
      name = "Mutation detections",
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Cell line dilution mutation detections",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      column_names_rot = 0,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          VAF_df[i, j], x, y, 
          gp = gpar(fontsize = 10, fontface = "bold", col = "#991425") )
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
  
  if (type == "dilution") {
    # add percentages to dilutions:
    colnames(supp_df)[!(colnames(supp_df) %in% c("Patient_id", "Site"))] <- paste0(
      colnames(supp_df)[!(colnames(supp_df) %in% c("Patient_id", "Site"))], "%" )
  }
  
  # make rownames patient ids and add pathology column:
  rownames(supp_df) <- supp_df$Patient_id
  supp_df$pathology <- " "
  supp_df <- subset(supp_df, select = -c(Patient_id, Site))
  
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
  
  # remake detection to include only fusions:
  supp_detect <- detect_df
  supp_detect[supp_df == " " & supp_detect == "detected"] <- "not detected"
  
  # create supporting reads heatmap:
  if (type == "patient") {
    sread_hm <- Heatmap(
      as.matrix(supp_detect), 
      name = "Fusion-supporting reads", 
      row_split = met_split,
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient EWSR1/FLI1 fusion-supporting reads",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          supp_df[i, j], x, y, 
          gp = gpar(fontsize = 12, fontface = "bold", col = "#430F82") )
      },
    )
  } else {
    sread_hm <- Heatmap(
      as.matrix(supp_detect), 
      name = "Fusion-supporting reads", 
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient EWSR1/FLI1 fusion-supporting reads",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      column_names_rot = 0,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          supp_df[i, j], x, y, 
          gp = gpar(fontsize = 12, fontface = "bold", col = "#430F82") )
      },
    )
  }

  # return heatmaps:
  return(list(VAF=VAF_hm, sread=sread_hm))

}
