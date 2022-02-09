fetch_sm_vafs <- function(sample_df, roi, gm_var, reverse_strand=FALSE ) {
  
  library(plyr)
  print(sample_df$Library_id)
  print(roi$name)
  
  # read in variants:
  try(
    vcf <- read.table(
      paste0(
        variant_path, "/", sample_df$Library_id, "/", 
        sample_df$Library_id, ".smCounter.anno.vcf" ),
      ), silent = TRUE)
  
  if (exists("vcf")) {
    
    # change TRUE to T where should be nucleotide:
    if (class(vcf$V4) == "logical") {
      vcf$V4 <- revalue(as.character(vcf$V4), c("TRUE" = "T")) 
    }
    
    # keep only those which passed filtering as gr, and add effect size and 
    # vaf columns:
    pass_var <- vcf[vcf$V7 == "PASS",]

    # convert to granges:
    pass_gr <- GRanges(
      seqnames = pass_var$V1,
      ranges = IRanges(start = pass_var$V2, end = pass_var$V2),
      ref = pass_var$V4,
      alt = pass_var$V5,
      qual = pass_var$V6,
      effect_size = sapply(strsplit(pass_var$V8, "\\|"), function(x) x[3]),
      UMT = gsub("^.*=", "", 
        gsub(",.*$", "", sapply(strsplit(pass_var$V8, ";"), function(x) x[4])) ),
      VMT = gsub("^.*=", "", 
        gsub(",.*$", "", sapply(strsplit(pass_var$V8, ";"), function(x) x[5])) ),
      VAF = round(as.numeric(
        gsub("^.*=", "", 
          gsub(",.*$", "", sapply(strsplit(pass_var$V8, ";"), function(x) x[6])) )
      )*100, 1),
      info = pass_var$V8
    )
    
    if (sample_df$Treatment == "tumour") {
      # remove tumour gDNA variants with VAF == 100 as likely germline:
      pass_gr <- pass_gr[pass_gr$VAF != 100]
    } else if (sample_df$Treatment != "100") {
      # remove ctDNA variants with VAF > 95 as likely germline:
      pass_gr <- pass_gr[pass_gr$VAF <= 95]
    }
    
    # convert effect size to numeric score:
    pass_gr$effect_size <- revalue(as.character(pass_gr$effect_size), 
      c("HIGH" = "3", "MODERATE" = "2", "LOW" = "1","MODIFIER" = "0") )

    # keep variants within roi:
    pint <- pintersect(pass_gr, roi)
    keep_var <- pint[pint$hit]
    
    # remove variants with low effect size:
    keep_var <- keep_var[keep_var$effect_size > 1]

    # remove benign germline variant chr17:7579472_G>C:
    top_var <- keep_var[
      !(start(keep_var) == 7579472 & keep_var$ref == "G" & keep_var$alt == "C") ]
    
#    # keep variant with highest effect size:
#    keep_var <- keep_var[which.max(keep_var$effect_size)]
#    top_var <- keep_var[which.max(keep_var$qual)]
    
    # change effect size back:
    top_var$effect_size <- revalue(as.character(top_var$effect_size), 
      c("3" = "HIGH", "2" = "MODERATE", "1" = "LOW", "0" = "MODIFIER" ) )

    # if reverse strand, convert nucleotides:
    if (reverse_strand) {
      top_var$ref <- chartr("ATGC","TACG", top_var$ref)
      top_var$alt <- chartr("ATGC","TACG", top_var$alt)
    }

    # remove germline variants:
    olaps <- findOverlaps(top_var, gm_var)
    for (i in seq_along(queryHits(olaps))) {
      if (top_var$ref[queryHits(olaps)[i]] == gm_var$ref[subjectHits(olaps)[i]] & 
        top_var$alt[queryHits(olaps)[i]] == gm_var$alt[subjectHits(olaps)[i]] ) {
        top_var <- top_var[-queryHits(olaps)[i]]
      }
    }
    
    # keep best quality variant call:
    if (length(top_var) > 0) {
      mcols(top_var) <- subset(mcols(top_var), select = -c(info, hit))
      return(as.data.frame(top_var))
    } else {
      return(data.frame(
        seqnames = "not_detected",
        start = "not_detected",
        end = "not_detected",
        width = "not_detected",
        strand = "not_detected",
        ref = "not_detected",
        alt = "not_detected",
        qual = "not_detected",
        effect_size = "not_detected",
        UMT = "not_detected",
        VMT = "not_detected",
        VAF = "not_detected" ))
    }
  } else {
    return(data.frame(
        seqnames = "not_detected",
        start = "not_detected",
        end = "not_detected",
        width = "not_detected",
        strand = "not_detected",
        ref = "not_detected",
        alt = "not_detected",
        qual = "not_detected",
        effect_size = "not_detected",
        UMT = "not_detected",
        VMT = "not_detected",
        VAF = "not_detected" ))
  }

}