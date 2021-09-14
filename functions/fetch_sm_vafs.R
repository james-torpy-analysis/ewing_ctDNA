fetch_sm_vafs <- function(sample_df, roi, alt_var) {
  
  print(sample_df)
  
  # read in variants:
  vcf <- read.table(
    paste0(
      variant_path, sample_df$Sample, "/", sample_df$Sample, 
      ".smCounter.anno.vcf")
  )
  
  # keep only those which passed filtering as gr, and add effect size and 
  # vaf columns:
  pass_var <- vcf[vcf$V7 == "PASS",]
  pass_gr <- GRanges(
    seqnames = pass_var$V1,
    ranges = IRanges(start = pass_var$V2, end = pass_var$V2),
    ref = pass_var$V4,
    alt = pass_var$V5,
    qual = pass_var$V6,
    effect_size = sapply(strsplit(pass_var$V8, "\\|"), function(x) x[3]),
    VAF = round(as.numeric(
      gsub("^.*=", "", sapply(strsplit(pass_var$V8, ";"), function(x) x[6]))
    )*100, 1),
    info = pass_var$V8
  )
  
  if (sample_df$Treatment == "tumour") {
    # remove tumour gDNA variants with VAF == 100 as likely germline:
    pass_gr <- pass_gr[pass_gr$VAF != 100]
  } else {
    # remove ctDNA variants with VAF > 95 as likely germline:
    pass_gr <- pass_gr[pass_gr$VAF <= 95]
  }
  
  # convert effect size to numeric score:
  pass_gr$effect_size[pass_gr$effect_size == "HIGH"] <- 3
  pass_gr$effect_size[pass_gr$effect_size == "MODERATE"] <- 2
  pass_gr$effect_size[pass_gr$effect_size == "LOW"] <- 1
  pass_gr$effect_size[pass_gr$effect_size == "MODIFIER"] <- 0
  
  # keep variant in roi with highest effect and top quality score, 
  # but remove those with 'MODIFIER' effect:
  top_var <- lapply(roi, function(x) {
    
    # keep variants within roi:
    pint <- pintersect(pass_gr, x)
    keep_var <- pint[pint$hit]
    
    # keep variants with highest effect size:
    keep_var <- keep_var[which.max(keep_var$effect_size)]
    keep_var[which.max(keep_var$qual)]
    
    # keep best quality variant call:
    if (length(keep_var) > 0) {
      if (keep_var$effect_size > 0) {
        return(keep_var)        
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
    
  })
  
  return(top_var)
}