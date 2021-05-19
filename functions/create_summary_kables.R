create_summary_kables <- function(results) {
  
  library(knitr)
  
  fusion_dfs <- list(
    patient = results$high_conf_bp$fusion_nos$patient_df,
    dilution = results$high_conf_bp$fusion_nos$dilution_df
  )
  
  # add less stringent calls to final df:
  fusion_dfs$patient$Less_stringent_fusion_detections <- 
    results$low_conf_bp$fusion_nos$patient_df$Detected_FLI1_EWSR1_fusion
  fusion_dfs$patient$Less_stringent_false_positives <- 
    results$low_conf_bp$fusion_nos$patient_df$False_EWSR1_fusions
  
  fusion_dfs$dilution$Less_stringent_fusion_detections <- 
    results$low_conf_bp$fusion_nos$dilution_df$Detected_FLI1_EWSR1_fusion
  fusion_dfs$dilution$Less_stringent_false_positives <- 
    results$low_conf_bp$fusion_nos$dilution_df$False_EWSR1_fusions
  
  # change high confidence colnames:
  colnames(fusion_dfs$patient) <- gsub(
    "Detected_FLI1_EWSR1_fusion", "Stringent_fusion_detections",
    colnames(fusion_dfs$patient)
  )
  colnames(fusion_dfs$patient) <- gsub(
    "False_EWSR1_fusions", "Stringent_false_positives",
    colnames(fusion_dfs$patient)
  )
  
  colnames(fusion_dfs$dilution) <- gsub(
    "Detected_FLI1_EWSR1_fusion", "Stringent_fusion_detections",
    colnames(fusion_dfs$dilution)
  )
  colnames(fusion_dfs$dilution) <- gsub(
    "False_EWSR1_fusions", "Stringent_false_positives",
    colnames(fusion_dfs$dilution)
  )
  
  # remove underscores and write to table:
  colnames(fusion_dfs$patient) <- gsub(
    "_", " ", 
    colnames(fusion_dfs$patient)
  )
  colnames(fusion_dfs$patient) <- gsub(
    "FLI1 EWSR1", "FLI1/EWSR1", 
    colnames(fusion_dfs$patient)
  )
  colnames(fusion_dfs$patient) <- gsub(
    "Treatment stage", "Treatment/stage", 
    colnames(fusion_dfs$patient)
  )
  
  colnames(fusion_dfs$dilution) <- gsub(
    "_", " ", 
    colnames(fusion_dfs$dilution)
  )
  
  return(
    list(
      patient = kable(
        fusion_dfs$patient, 
        align = rep("c", ncol(fusion_dfs$patient))
      ),
      dilution = kable(
        fusion_dfs$dilution, 
        align = rep("c", ncol(fusion_dfs$dilution))
      )
    )
  )
  
}