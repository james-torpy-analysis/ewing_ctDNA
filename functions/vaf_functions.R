reverseStrand <- function (x)
{

    old2new <- c("+" = "-", "-" = "+", "*" = "*")
    strand(x) <- old2new[as.character(strand(x))]
    return(x)

}

writeSam <- function (bam_in, selected, sam_out)
{

    ## read bam alignments
    ga <- readGAlignments(bam_in, use.names = TRUE, param = param)
    cols <- c(
        "qname",
        "flag",
        "rname",
        "pos",
        "mapq",
        "cigar",
        "mrnm",
        "mpos",
        "isize",
        "seq",
        "qual")

    ## subset to select reads
    index <- which(mcols(ga)$qname %in% selected)
    sam <- mcols(ga)[index, cols]

    ## write sam
    if (file.exists(sam_out)) file.remove(sam_out)

    system(paste0("samtools view -H ", bam_in, " > ", sam_out))
    write.table(sam,
                file = sam_out,
                append = TRUE,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    
    # convert to bam and index:
    bam_out <- gsub("sam", "bam", sam_out)
    system(paste0("samtools view -bh ", sam_out, " > ", bam_out))
    system(paste0("samtools index ", bam_out))
    
    # remove sams:
    system(paste0("rm ", sam_out, "*"))

}

filter_and_VAF <- function(fusions_gr, gene_a_name, gene_b_name, hist_dir, 
  out_bam_dir ) {

  for (i in seq_along(fusions_gr)) {
    fusion <- fusions_gr[i]

    gene_a_breakpoint <- fusion    
    gene_b_breakpoint <- GRanges(fusion$join_chr,
      IRanges(fusion$join_coord, fusion$join_coord),fusion$join_strand )
    mcols(gene_a_breakpoint) <- NULL
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
    supp_fwd_rev <- intersect(
      names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
      names(which(sum(width(pintersect(R2, gene_b_dnstream_rev))) >= min_overlap_R2)))
    supp_fwd_fwd <- intersect(
      names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
      names(which(sum(width(pintersect(R2, gene_b_upstream))) >= min_overlap_R2)))

    ## diagnost plots for overlap
    R1_nonsupp_fwd <- R1[which(names(R1) %in% nonsupp_fwd)]
    R2_nonsupp_fwd <- R2[which(names(R2) %in% nonsupp_fwd)]
    R1_supp_fwd_rev <- R1[which(names(R1) %in% supp_fwd_rev)]
    R2_supp_fwd_rev <- R2[which(names(R2) %in% supp_fwd_rev)]

    pdf(file.path(hist_dir, paste0("hist_overlap_", breakpoint, "_fwd.pdf")))
    par(mfrow = c(2, 2))
    if (length(pintersect(R1_nonsupp_fwd, gene_a_upstream))>0) {
      hist(-sum(width(pintersect(R1_nonsupp_fwd, gene_a_upstream))),
        xlim = c(-180, 0), xlab = "Overlap [bp]",
        main = paste0("fusion non-supporting ", gene_a_name, " upstream") )}
    if (length(pintersect(R2_nonsupp_fwd, gene_a_dnstream_rev))>0) {
      hist(sum(width(pintersect(R2_nonsupp_fwd, gene_a_dnstream_rev))),
         xlim = c(0, 180), xlab = "Overlap [bp]",
         main = paste0("fusion non-supporting ", gene_a_name, " dnstream") )}
    if (length(pintersect(R1_supp_fwd_rev, gene_a_upstream))>0) {
      hist(-sum(width(pintersect(R1_supp_fwd_rev, gene_a_upstream))),
           xlim = c(-180, 0), xlab = "Overlap [bp]",
           main = paste0("fusion supporting ", gene_a_name, " upstream") )}
    if (length(pintersect(R1_supp_fwd_rev, gene_a_upstream))>0) {
      hist(sum(width(pintersect(R1_supp_fwd_rev, gene_a_upstream))),
           xlim = c(0, 180), xlab = "Overlap [bp]",
           main = paste0("fusion supporting ", gene_b_name, " dnstream") )}
    dev.off()

    ## non-supporting reads that satisfy overlap criteria for reciprocal fusion
    nonsupp_rev <- intersect(
      names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
      names(which(sum(width(pintersect(R2, gene_a_upstream))) >= min_overlap_R2)))
    ## supporting reads that satisfy overlap criteria for reciprocal fusion
    supp_rev_fwd <- intersect(
      names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
      names(which(sum(width(pintersect(R2, gene_b_upstream))) >= min_overlap_R2)))
    supp_rev_rev <- intersect(
      names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
      names(which(sum(width(pintersect(R2, gene_b_dnstream_rev))) >= min_overlap_R2)))
    
    ## diagnost plots for overlap
    R1_nonsupp_rev <- R1[which(names(R1) %in% nonsupp_rev)]
    R2_nonsupp_rev <- R2[which(names(R2) %in% nonsupp_rev)]
    R1_supp_rev <- R1[which(names(R1) %in% supp_rev)]
    R2_supp_rev <- R2[which(names(R2) %in% supp_rev)]
    
    pdf(file.path(hist_dir, paste0("hist_overlap_", breakpoint, "_rev.pdf")))
    par(mfrow = c(2, 2))
    if (length(pintersect(R2_nonsupp_rev, gene_a_upstream))>0) {
      hist(-sum(width(pintersect(R2_nonsupp_rev, gene_a_upstream))),
        xlim = c(-180, 0), xlab = "Overlap [bp]",
        main = paste0("reciproc non-supporting ", gene_a_name, " upstream") )}
    if (length(pintersect(R1_nonsupp_rev, gene_a_dnstream_rev))>0) {
      hist(sum(width(pintersect(R1_nonsupp_rev, gene_a_dnstream_rev))),
        xlim = c(0, 180), xlab = "Overlap [bp]",
        main = paste0("reciproc non-supporting ", gene_a_name, " dnstream") )}
    if (length(pintersect(R2_supp_rev, gene_b_upstream))>0) {
      hist(-sum(width(pintersect(R2_supp_rev, gene_b_upstream))),
          xlim = c(-180, 0), xlab = "Overlap [bp]",
          main = paste0("reciproc supporting ", gene_b_name, " upstream") )}
    if (length(pintersect(R1_supp_rev, gene_a_dnstream_rev))>0) {
      hist(sum(width(pintersect(R1_supp_rev, gene_a_dnstream_rev))),
          xlim = c(0, 180), xlab = "Overlap [bp]",
          main = paste0("reciproc supporting ", gene_a_name, " dnstream") )}
    dev.off()
  
    ## calculate VAFs for all translocations
    print(length(nonsupp_fwd))
    print(length(supp_fwd_rev))
    VAF_fwd_rev <- length(supp_fwd_rev) / (length(nonsupp_fwd) + length(supp_fwd_rev))
    print(VAF_fwd_rev)

    print(length(nonsupp_fwd))
    print(length(supp_fwd_fwd))
    VAF_fwd_fwd <- length(supp_fwd_fwd) / (length(nonsupp_fwd) + length(supp_fwd_fwd))
    print(VAF_fwd_fwd)
    
    print(length(nonsupp_rev))
    print(length(supp_rev_fwd))
    VAF_rev_fwd <- length(supp_rev_fwd) / (length(nonsupp_rev) + length(supp_rev_fwd))
    print(VAF_rev_fwd)

    print(length(nonsupp_rev))
    print(length(supp_rev_rev))
    VAF_rev_rev <- length(supp_rev_rev) / (length(nonsupp_rev) + length(supp_rev_rev))
    print(VAF_rev_rev)
    
    
    if (i==1) {
      VAFs <- list(data.frame(VAF_fwd_rev, VAF_fwd_fwd, VAF_rev_fwd, VAF_rev_rev))
    } else {
      VAFs[[i]] <- data.frame(VAF_fwd_rev, VAF_fwd_fwd, VAF_rev_fwd, VAF_rev_rev) }
    names(VAFs)[i] <- paste0("fusion_", start(fusions_gr[i]))

    VAFs[[i]]$Fusion <- paste0(
      gsub(":\\+", "-", as.character(fusion)), fusion$join_chr, ":", 
      fusion$join_coord )
    VAFs[[i]]$Forward_non_supporting <- length(nonsupp_fwd)
    VAFs[[i]]$Forward_reverse_supporting <- length(supp_fwd_rev)
    VAFs[[i]]$Forward_reverse_total <- VAFs[[i]]$Forward_non_supporting + 
      VAFs[[i]]$Forward_reverse_supporting
    VAFs[[i]]$Forward_forward_supporting <- length(supp_fwd_fwd)
    VAFs[[i]]$Forward_forward_total <- VAFs[[i]]$Forward_non_supporting + 
      VAFs[[i]]$Forward_forward_supporting
    VAFs[[i]]$Reverse_non_supporting <- length(nonsupp_rev)
    VAFs[[i]]$Reverse_forward_supporting <- length(supp_rev_fwd)
    VAFs[[i]]$Reverse_forward_total <- VAFs[[i]]$Reverse_non_supporting + 
      VAFs[[i]]$Reverse_forward_supporting
    VAFs[[i]]$Reverse_reverse_supporting <- length(supp_rev_rev)
    VAFs[[i]]$Reverse_reverse_total <- VAFs[[i]]$Reverse_non_supporting + 
      VAFs[[i]]$Reverse_reverse_supporting
    VAFs[[i]]$All_supporting <- VAFs[[i]]$Forward_reverse_supporting + 
      VAFs[[i]]$Forward_forward_supporting + VAFs[[i]]$Reverse_forward_supporting + 
      VAFs[[i]]$Reverse_reverse_supporting

    ## write reads to SAM for inspection
    writeSam(file_bam, nonsupp_fwd, file.path(out_bam_dir, paste0("reads_", breakpoint, "_nonsupp_fwd.sam")))
    writeSam(file_bam, supp_fwd_rev, file.path(out_bam_dir, paste0("reads_", breakpoint, "_supp_fwd_rev.sam")))
    writeSam(file_bam, supp_fwd_fwd, file.path(out_bam_dir, paste0("reads_", breakpoint, "_supp_fwd_fwd.sam")))
    writeSam(file_bam, nonsupp_rev, file.path(out_bam_dir, paste0("reads_", breakpoint, "_nonsupp_rev.sam")))
    writeSam(file_bam, supp_rev_fwd, file.path(out_bam_dir, paste0("reads_", breakpoint, "_supp_rev_fwd.sam")))  
    writeSam(file_bam, supp_rev_rev, file.path(out_bam_dir, paste0("reads_", breakpoint, "_supp_rev_rev.sam")))  
  }
  return(VAFs)
}

