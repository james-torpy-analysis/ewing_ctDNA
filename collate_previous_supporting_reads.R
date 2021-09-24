args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]

#projectname <- "ewing_ctDNA"
#samplename <- "409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
result_dir <- paste0(project_dir, "results_200821/")
in_path <- paste0(result_dir, "VAF_calculation/", samplename, "/Rdata/")

load(file.path(in_path, "VAFs_calculated.Rdata"))

result_dir <- paste0(project_dir, "results_200821/")
out_path <- paste0(result_dir, "VAF_calculation/", samplename, "/")

for (i in seq_along(fusions)) {
  
  fusion <- fusions[i]
  
  gene_a_breakpoint <- GRanges(
    fusion$join_chr,
    IRanges(
      fusion$join_coord,
      fusion$join_coord),
    fusion$join_strand)
  gene_b_breakpoint <- fusion
  mcols(gene_b_breakpoint) <- NULL
  
  breakpoint <- gsub(":", "_", as.character(gene_a_breakpoint), fixed = TRUE)
  print(breakpoint)
  
  ## consider 1Mb upstream or downstream of EWSR1 breakpoint
  gene_a_upstream <- flank(gene_a_breakpoint, 1e6, start = TRUE)
  gene_a_upstream <- resize(gene_a_upstream, 1e6 + 1) ## add breakpoint position
  gene_a_dnstream <- flank(gene_a_breakpoint, 1e6, start = FALSE)
  gene_a_dnstream_rev <- reverseStrand(gene_a_dnstream)
  
  ## consider 1Mb upstream or downstream of fusion partner breakpoint
  gene_b_upstream <- flank(gene_b_breakpoint, 1e6, start = TRUE)
  gene_b_upstream <- resize(gene_b_upstream, 1e6 + 1) ## add breakpoint position
  gene_b_dnstream <- flank(gene_b_breakpoint, 1e6, start = FALSE)
  gene_b_dnstream_rev <- reverseStrand(gene_b_dnstream)
  
  ## non-supporting reads that satisfy overlap criteria for driver fusion
  nonsupp_fwd <- intersect(
    names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
    names(which(sum(width(pintersect(R2, gene_a_dnstream_rev))) >= min_overlap_R2)))
  ## supporting reads that satisfy overlap criteria for driver fusion
  supp_fwd <- intersect(
    names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
    names(which(sum(width(pintersect(R2, gene_b_dnstream_rev))) >= min_overlap_R2)))
  
  ## diagnost plots for overlap
  
  R1_nonsupp_fwd <- R1[which(names(R1) %in% nonsupp_fwd)]
  R2_nonsupp_fwd <- R2[which(names(R2) %in% nonsupp_fwd)]
  R1_supp_fwd <- R1[which(names(R1) %in% supp_fwd)]
  R2_supp_fwd <- R2[which(names(R2) %in% supp_fwd)]
  
  ## non-supporting reads that satisfy overlap criteria for reciprocal fusion
  nonsupp_rev <- intersect(
    names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
    names(which(sum(width(pintersect(R2, gene_a_upstream))) >= min_overlap_R2)))
  ## supporting reads that satisfy overlap criteria for reciprocal fusion
  supp_rev <- intersect(
    names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
    names(which(sum(width(pintersect(R2, gene_b_upstream))) >= min_overlap_R2)))
  
  ## diagnost plots for overlap
  
  R1_nonsupp_rev <- R1[which(names(R1) %in% nonsupp_rev)]
  R2_nonsupp_rev <- R2[which(names(R2) %in% nonsupp_rev)]
  R1_supp_rev <- R1[which(names(R1) %in% supp_rev)]
  R2_supp_rev <- R2[which(names(R2) %in% supp_rev)]
  
  ## calculate VAFs for both translocations
  print(length(nonsupp_fwd))
  print(length(supp_fwd))
  VAF_fwd <- length(supp_fwd) / (length(nonsupp_fwd) + length(supp_fwd))
  print(VAF_fwd)
  
  print(length(nonsupp_rev))
  print(length(supp_rev))
  VAF_rev <- length(supp_rev) / (length(nonsupp_rev) + length(supp_rev))
  print(VAF_rev)
  
  if (i==1) {
    VAFs <- list(
      data.frame(VAF_fwd, VAF_rev, length(supp_fwd), length(supp_rev))
    )
    VAFs[[i]]$supp_total <- VAFs[[i]]$length.supp_fwd. + VAFs[[i]]$length.supp_rev.
  } else {
    VAFs[[i]] <- data.frame(VAF_fwd, VAF_rev, length(supp_fwd), length(supp_rev))
    VAFs[[i]]$supp_total <- VAFs[[i]]$length.supp_fwd. + VAFs[[i]]$length.supp_rev.
  }
  names(VAFs)[i] <- paste0("fusion_", fusions[i]$join_coord)
}

VAF_df <- do.call("rbind", VAFs)
  VAF_df$conf <- fusions$conf
  
write.table(
  VAF_df,
  paste0(out_path, "VAFs_with_supp_read_no.tsv"),
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T)
