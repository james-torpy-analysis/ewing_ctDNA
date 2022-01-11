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

    ## remove reads without mapped mates:
    sam$mpos[is.na(sam$mpos)] <- "0"
    sam$isize[is.na(sam$isize)] <- "0"

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

filter_and_VAF <- function(fusions_gr, max_suppl_dist, min_suppl_reads,
  gene_a_name, gene_b_name, R1, R2, hist_dir, out_bam_dir ) {

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

    ## non-supporting reads
    nonsupp <- list(
      up2dn=intersect(
        names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_a_dnstream_rev))) >= min_overlap_R2)) ),
      dn2up=intersect(
        names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_a_upstream))) >= min_overlap_R2)) ))

    ## supporting reads
    supp <- list(
      up2dn=intersect(
        names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_b_dnstream_rev))) >= min_overlap_R2)) ),
      up2up=intersect(
        names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_b_upstream))) >= min_overlap_R2)) ),
      dn2up=intersect(
        names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_b_upstream))) >= min_overlap_R2)) ),
      dn2dn=intersect(
        names(which(sum(width(pintersect(R1, gene_a_dnstream_rev))) >= min_overlap_R1)),
        names(which(sum(width(pintersect(R2, gene_b_dnstream_rev))) >= min_overlap_R2)) ))

    # ensure supplementary reads fall within n bp of breakpoint:
    fusion2 <- GRanges(
      seqnames=fusion$join_chr,
      ranges=IRanges(start=fusion$join_coord, width=1),
      strand=fusion$join_strand )

    supp <- lapply(supp, function(x) {
      # fetch all mates/supplementary aligns of each supporting read:
      reads <- c(unlist(R1[names(R1) %in% x]), unlist(R2[names(R2) %in% x]))
      # determine distance between reads and closest breakpoint:
      strand(reads) <- "*"
      reads$dist_from_bp <- distance(reads, fusion)
      reads$dist_from_bp2 <- distance(reads, fusion2)
      reads$dist_from_bp[is.na(reads$dist_from_bp)] <- reads$dist_from_bp2[
        is.na(reads$dist_from_bp) ]
      # remove supplementary (split) reads < max_suppl_dist from breakpoint:
      remove <- unique(gsub("\\..*$", "", names(reads)[
        reads$flag >=2048 & reads$dist_from_bp > max_suppl_dist ]))
      if (length(remove)>0) {
        return(x[x!=remove])
      } else {
        return(x)
      } 
    })

    # require either 1 discordant read pair or 2 supplementary reads to keep:
    supp <- lapply(supp, function(x) {
      if (length(x)==1) {
        reads <- c(unlist(R1[names(R1) %in% x]), unlist(R2[names(R2) %in% x]))
        if (length(unique(seqnames(reads[reads$flag <2048]))) <2 & 
          length(which(reads$flag >=2048)) <min_suppl_reads) {
          return(as.character(NULL))
        } else {return(x)}
      } else {return(x)} })
      

    ## diagnost plots for overlap
    R1_nonsupp_up2dn <- R1[which(names(R1) %in% nonsupp$up2dn)]
    R2_nonsupp_up2dn <- R2[which(names(R2) %in% nonsupp$up2dn)]
    R1_supp_up2dn <- R1[which(names(R1) %in% supp$up2dn)]
    R2_supp_up2dn <- R2[which(names(R2) %in% supp$up2dn)]
    R1_nonsupp_dn2up <- R1[which(names(R1) %in% nonsupp$dn2up)]
    R2_nonsupp_dn2up <- R2[which(names(R2) %in% nonsupp$dn2up)]
    R1_supp_dn2up <- R1[which(names(R1) %in% supp$dn2up)]
    R2_supp_dn2up <- R2[which(names(R2) %in% supp$dn2up)]

    pdf(file.path(hist_dir, paste0("hist_overlap_", breakpoint, "_up2dn.pdf")))
    par(mfrow = c(2, 2))
    if (length(pintersect(R1_nonsupp_up2dn, gene_a_upstream))>0) {
      hist(-sum(width(pintersect(R1_nonsupp_up2dn, gene_a_upstream))),
        xlim = c(-180, 0), xlab = "Overlap [bp]",
        main = paste0("fusion non-supporting ", gene_a_name, " upstream") )}
    if (length(pintersect(R2_nonsupp_up2dn, gene_a_dnstream_rev))>0) {
      hist(sum(width(pintersect(R2_nonsupp_up2dn, gene_a_dnstream_rev))),
         xlim = c(0, 180), xlab = "Overlap [bp]",
         main = paste0("fusion non-supporting ", gene_a_name, " dnstream") )}
    if (length(pintersect(R1_supp_up2dn, gene_a_upstream))>0) {
      hist(-sum(width(pintersect(R1_supp_up2dn, gene_a_upstream))),
           xlim = c(-180, 0), xlab = "Overlap [bp]",
           main = paste0("fusion supporting ", gene_a_name, " upstream") )}
    if (length(pintersect(R1_supp_up2dn, gene_a_upstream))>0) {
      hist(sum(width(pintersect(R1_supp_up2dn, gene_a_upstream))),
           xlim = c(0, 180), xlab = "Overlap [bp]",
           main = paste0("fusion supporting ", gene_b_name, " dnstream") )}
    dev.off()
    pdf(file.path(hist_dir, paste0("hist_overlap_", breakpoint, "_dn2up.pdf")))
    par(mfrow = c(2, 2))
    if (length(pintersect(R2_nonsupp_dn2up, gene_a_upstream))>0) {
      hist(-sum(width(pintersect(R2_nonsupp_dn2up, gene_a_upstream))),
        xlim = c(-180, 0), xlab = "Overlap [bp]",
        main = paste0("reciproc non-supporting ", gene_a_name, " upstream") )}
    if (length(pintersect(R1_nonsupp_dn2up, gene_a_dnstream_rev))>0) {
      hist(sum(width(pintersect(R1_nonsupp_dn2up, gene_a_dnstream_rev))),
        xlim = c(0, 180), xlab = "Overlap [bp]",
        main = paste0("reciproc non-supporting ", gene_a_name, " dnstream") )}
    if (length(pintersect(R2_supp_dn2up, gene_b_upstream))>0) {
      hist(-sum(width(pintersect(R2_supp_dn2up, gene_b_upstream))),
          xlim = c(-180, 0), xlab = "Overlap [bp]",
          main = paste0("reciproc supporting ", gene_b_name, " upstream") )}
    if (length(pintersect(R1_supp_dn2up, gene_a_dnstream_rev))>0) {
      hist(sum(width(pintersect(R1_supp_dn2up, gene_a_dnstream_rev))),
          xlim = c(0, 180), xlab = "Overlap [bp]",
          main = paste0("reciproc supporting ", gene_a_name, " dnstream") )}
    dev.off()

    ## calculate VAFs for all translocations
    print(length(nonsupp$up2dn))
    print(length(supp$up2dn))
    VAF_up2dn <- length(supp$up2dn) / (length(nonsupp$up2dn) + length(supp$up2dn))
    print(VAF_up2dn)

    print(length(nonsupp$up2dn))
    print(length(supp$up2up))
    VAF_up2up <- length(supp$up2up) / (length(nonsupp$up2dn) + length(supp$up2up))
    print(VAF_up2up)
    
    print(length(nonsupp$dn2up))
    print(length(supp$dn2up))
    VAF_dn2up <- length(supp$dn2up) / (length(nonsupp$dn2up) + length(supp$dn2up))
    print(VAF_dn2up)

    print(length(nonsupp$dn2up))
    print(length(supp$dn2dn))
    VAF_dn2dn <- length(supp$dn2dn) / (length(nonsupp$dn2up) + length(supp$dn2dn))
    print(VAF_dn2dn)
    
    
    if (i==1) {
      VAFs <- list(round(data.frame(VAF_forward=VAF_up2dn, VAF_up_to_upstream=VAF_up2up, 
        VAF_reverse=VAF_dn2up, VAF_down_to_downstream=VAF_dn2dn ), 4))
    } else {
      VAFs[[i]] <- round(data.frame(VAF_forward=VAF_up2dn, VAF_up_to_upstream=VAF_up2up, 
        VAF_reverse=VAF_dn2up, VAF_down_to_downstream=VAF_dn2dn ), 4) }
    names(VAFs)[i] <- paste0("fusion_", start(fusions_gr[i]))

    VAFs[[i]]$Fusion <- paste0(gsub(":\\+", "-", as.character(fusion)), 
      fusion$join_chr, ":", fusion$join_coord )
    VAFs[[i]]$Forward_non_supporting <- length(nonsupp$up2dn)
    VAFs[[i]]$Forward_supporting <- length(supp$up2dn)
    VAFs[[i]]$Forward_total <- VAFs[[i]]$Forward_non_supporting + 
      VAFs[[i]]$Forward_supporting
    VAFs[[i]]$Up_to_upstream_supporting <- length(supp$up2up)
    VAFs[[i]]$Up_to_upstream_total <- VAFs[[i]]$Forward_non_supporting + 
      VAFs[[i]]$Up_to_upstream_supporting
    VAFs[[i]]$Reverse_non_supporting <- length(nonsupp$dn2up)
    VAFs[[i]]$Down_to_upstream_supporting <- length(supp$dn2up)
    VAFs[[i]]$Down_to_upstream_total <- VAFs[[i]]$Reverse_non_supporting + 
      VAFs[[i]]$Down_to_upstream_supporting
    VAFs[[i]]$Down_to_downstream_supporting <- length(supp$dn2dn)
    VAFs[[i]]$Down_to_downstream_total <- VAFs[[i]]$Reverse_non_supporting + 
      VAFs[[i]]$Down_to_downstream_supporting
    VAFs[[i]]$All_supporting <- VAFs[[i]]$Forward_supporting + 
      VAFs[[i]]$Up_to_upstream_supporting + VAFs[[i]]$Down_to_upstream_supporting + 
      VAFs[[i]]$Down_to_downstream_supporting

    ## write reads to SAM for inspection
    writeSam(file_bam, nonsupp$up2dn, file.path(out_bam_dir, 
      paste0("reads_", breakpoint, "_nonsupp_fwd.sam") ))
    writeSam(file_bam, supp$up2dn, file.path(out_bam_dir, 
      paste0("reads_", breakpoint, "_supp_fwd.sam") ))
    writeSam(file_bam, supp$up2up, file.path(out_bam_dir, 
      paste0("reads_", breakpoint, "_supp_up2up.sam") ))
    writeSam(file_bam, nonsupp$dn2up, file.path(out_bam_dir, 
      paste0("reads_", breakpoint, "_nonsupp_rev.sam") ))
    writeSam(file_bam, supp$dn2up, file.path(out_bam_dir, 
      paste0("reads_", breakpoint, "_supp_rev.sam") ))  
    writeSam(file_bam, supp$dn2dn, file.path(out_bam_dir, 
      paste0("reads_", breakpoint, "_supp_dn2dn.sam") ))
  }
  return(VAFs)
}

