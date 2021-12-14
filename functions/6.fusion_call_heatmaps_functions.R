mutation_heatmap <- function(
  fusion_df,
  type,
  detection_cols,
  path_cols
) {
  
  library(naturalsort)
  library(reshape2)
  library(ComplexHeatmap)
  library(tibble)
  library(dplyr)
  
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
    x <- x[grep("VAF", names(x))]
    x <- x[grep("GeneGlobe", names(x), invert=T)]
    x <- x[grep("ddPCR", names(x), invert=T)]
    x <- x[grep("not_detected", x, invert=T)]
    
    if (length(x) > 0) {
      # round VAFs:
      VAF <- round2(as.numeric(x), 1)
      names(VAF) <- names(x)
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
  VAF_df$mutation_VAF <- gsub("_forward|_reverse", "", VAF_df$mutation_VAF)
  
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
  VAF_df$fusion_pathology <- ""
  VAF_df$pointmut_pathology <- ""
  VAF_df <- VAF_df[,colnames(VAF_df) %in% 
      c("fusion_pathology", "pointmut_pathology", 
      unique(as.character(fusion_df$Treatment.dilution)) )]
  VAF_df <- VAF_df %>%
    select(fusion_pathology, pointmut_pathology, everything())
  
  # merge pathology columns:
  fusion_path <- fusion_df[, grep("pathology", colnames(fusion_df))]
  fusion_df$fusion_pathology <- apply(fusion_path, 1, function(x) {
    fus <- names(x)[which(x != "not_detected")]
    if (length(fus) == 1) {return(fus)} else {return(x[1])} })
  fusion_df$fusion_pathology <- gsub("rearrangement|translocation", 
    "SV_pathology", fusion_df$fusion_pathology)
  
  pointmut_path <- fusion_df[, grep("Sanger", colnames(fusion_df))]
  fusion_df$pointmut_pathology <- apply(pointmut_path, 1, function(x) {
    pointmut <- names(x)[which(x != "not_detected")]
    if (length(pointmut) == 1) {
      return(paste0(gsub("Sanger_|_point_mut", "", pointmut), "_pathology"))
    } else if (length(pointmut) == 0) {return(x[1])} else {
      return("both_point_mut_pathology")
    } })

  # merge pathology results with dfs:
  path_df <- subset(fusion_df, select = c(Patient_id, fusion_pathology, 
    pointmut_pathology ))
  path_df <- path_df[!duplicated(path_df),]

  # bind to detect_df:
  detect_df <- merge(path_df, detect_df, by="Patient_id")
  detect_df <- column_to_rownames(detect_df, "Patient_id")

  # order dfs by name, pathology status, and site:
  if (type == "patient") {
    detect_df <- as.data.frame(
      detect_df[naturalorder(gsub("^.*_|M0|P", "", rownames(detect_df))),] )
    detect_df <- arrange(detect_df, fusion_pathology, pointmut_pathology)
    detect_df$Site <- factor(detect_df$Site, levels=c("primary", "met", "unknown"))
    detect_df <- detect_df[order(detect_df$Site),]
  } else {detect_df <- detect_df[order(rownames(detect_df)),]}
  
  # define site(row) and pathology(column) split vectors:
  if (type == "patient") {met_split <- detect_df$Site}
  path_split <- factor(
    c(rep("path", 2), rep("non_path", length(unique(fusion_df$Treatment.dilution)))),
    levels = c("path", "non_path") )
  detect_df <- subset(detect_df, select = -Site)
  
  VAF_df <- VAF_df[match(rownames(detect_df), rownames(VAF_df)),]
  
  # make NAs 'no_sample':
  detect_df <- as.data.frame(apply(detect_df, 2, function(x) {
    x[is.na(x)] <- "no_sample"
    return(x)
  }))
  
  # name detection colours:
  path_names <- unique(c(detect_df$fusion_pathology, detect_df$pointmut_pathology))
  path_names <- path_names[grep("pathology", path_names)]
  path_cols <- path_cols[1:length(path_names)]
  names(path_cols) <- path_names
  hm_cols <- c(detection_cols, path_cols)
  
  # remove underscores:
  detect_df <- as.data.frame(gsub("_", " ", as.matrix(detect_df)))
  rownames(detect_df) <- gsub("_", " ", rownames(detect_df))
  colnames(detect_df) <- gsub("_", " ", colnames(detect_df))
  names(hm_cols) <- gsub("_", " ", names(hm_cols))
  
  # make NAs blank:
  VAF_df <- apply(VAF_df, 2, function(x) {
    x[is.na(x)] <- ""
    return(x) })

  if (type == "dilution") {
    # order by dilution level:
    detect_df <- detect_df[, 
      c("fusion pathology", "pointmut pathology", unique(fusion_df$Treatment.dilution)) ]
    VAF_df <- VAF_df[,
      c("fusion_pathology", "pointmut_pathology", unique(fusion_df$Treatment.dilution)) ]
    # add percentages to dilutions:
    colnames(detect_df)[colnames(detect_df) != "fusion pathology" & 
      colnames(detect_df) != "pointmut pathology"] <- paste0(
      colnames(detect_df)[colnames(detect_df) != "fusion pathology" & 
      colnames(detect_df) != "pointmut pathology"], "%" )
    colnames(VAF_df)[colnames(VAF_df) != "fusion pathology" & 
      colnames(VAF_df) != "pointmut pathology"] <- paste0(
      colnames(VAF_df)[colnames(VAF_df) != "fusion pathology" & 
      colnames(VAF_df) != "pointmut pathology"], "%" ) }
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
      }
    )
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
  
  # merge pathology columns:
  temp_supp <- fusion_df[, grep("forward_supporting", colnames(fusion_df))]
  fusion_df$all_forward_supporting <- apply(temp_supp, 1, function(x) {
    supp <- x[which(x != "not_detected")]
    if (length(supp) == 1) {
      return(paste0(gsub("_forward_supporting", "", names(supp)), ": ", supp))
    } else if (length(supp) == 0) {return(x[1])
    } else {
      print("Error: multiple fusions detected when there should only be one, check data")
    } })

  # create supporting read annotation:
  supp_df <- subset(fusion_df, select = c(Patient_id, Site, Treatment.dilution,
    fusion_pathology, all_forward_supporting ))
  
  # make NA or not detected values blank:
  supp_df$all_forward_supporting[
    is.na(supp_df$all_forward_supporting) | 
    supp_df$all_forward_supporting == "not_detected" ] <- " "
  
  # cast to wide format:
  supp_df <- dcast(supp_df, Patient_id + Site + fusion_pathology ~Treatment.dilution, 
    value.var = "all_forward_supporting" )
  
  if (type == "dilution") {
    # add percentages to dilutions:
    colnames(supp_df)[!(colnames(supp_df) %in% 
      c("Patient_id", "Site", "fusion_pathology"))] <- paste0(
      colnames(supp_df)[!(colnames(supp_df) %in% 
        c("Patient_id", "Site", "fusion_pathology"))], "%" ) }
  
  # make rownames patient ids and add pathology column:
  rownames(supp_df) <- supp_df$Patient_id
  supp_df <- subset(supp_df, select = -c(Patient_id, Site))
  
  # remove underscores:
  rownames(supp_df) <- gsub("_", " ", rownames(supp_df))
  colnames(supp_df) <- gsub("_", " ", colnames(supp_df))
  
  # order df:
  supp_df <- supp_df[match(rownames(detect_df), rownames(supp_df)), 
    match(colnames(detect_df)[
      grep("pointmut pathology", colnames(detect_df), invert=T )], colnames(supp_df)) ]
  
  # make NAs and pathology columns blank:
  supp_df <- apply(supp_df, 2, function(x) {
    x[is.na(x)] <- " "
    x[grep("pathology", x)] <- " "
    return(x) })

  # remake detection to include only fusions:
  supp_detect <- detect_df[,grep("pointmut pathology", colnames(detect_df), invert=T)]
  supp_detect[supp_df == " " & supp_detect == "detected"] <- "not detected"
  
  # create supporting reads heatmap:
  if (type == "patient") {
    sread_hm <- Heatmap(
      as.matrix(supp_detect), 
      name = "Fusion-supporting reads", 
      row_split = met_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient fusion-supporting reads",
      column_title_gp = gpar(fontsize = 16, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          supp_df[i, j], x, y, 
          gp = gpar(fontsize = 12, fontface = "bold", col = "#430F82") )})
  } else {
    sread_hm <- Heatmap(
      as.matrix(supp_detect), 
      name = "Fusion-supporting reads",
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Cell line fusion-supporting reads",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      column_names_rot = 0,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          supp_df[i, j], x, y, 
          gp = gpar(fontsize = 12, fontface = "bold", col = "#430F82") )})}


  # return heatmaps:
  return(list(VAF=VAF_hm, sread=sread_hm))

}
