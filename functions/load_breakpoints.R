load_breakpoints <- function(samplename, in_dir, high_conf = TRUE) {

  # load VCF and convert to granges:
  print(paste0("Loading ", samplename, " VCF..."))

  if (high_conf) {
    in_file <- paste0(in_dir, samplename, ".svaba.sv.vcf")
  } else {
    in_file <- paste0(in_dir, samplename, ".svaba.semifiltered.sv.formatted.vcf")
  }

  vcf_df <- tryCatch(
    read.table(
      in_file,
      sep = "\t"
    ),
    warning = function(w) {NA},
    error = function(e) {NA}
  )
  
  if (any(!is.na(vcf_df))) {

    # remove alternate chromosomes:
    vcf_df <- vcf_df[grep("chrG|chrK|chrJ", vcf_df$V1, invert = T),]
    vcf_df <- vcf_df[grep("chrG|chrK|chrJ", vcf_df$V5, invert = T),]

    if (is.data.frame(vcf_df)) {
      
      # fetch all terms:
      info_terms <- unique(gsub("=.*$", "", unlist(strsplit(vcf_df$V8, ";"))))
      
      term_df <- data.frame(
        info_terms = info_terms
      )
      
      # split additional info field:
      for (i in 1:nrow(vcf_df)) {
        
        # split categories:
        spl <- strsplit(vcf_df$V8[i], ";")[[1]]
        
        # make IMPRECISE values true if they exist:
        spl[spl=="IMPRECISE"] <- "IMPRECISE=TRUE"
        
        # split values from categories:
        spl2 <- strsplit(spl, "=")
        temp_df <- as.data.frame(do.call("rbind", spl2))
        colnames(temp_df) <- c("info_terms", i)
        
        # merge all terms df with temp_df:
        info_df <- merge(term_df, temp_df, by = "info_terms")
        
        # merge with additional_info_df:
        if (i==1) {
          additional_info <- info_df
        } else {
          additional_info <- merge(
            additional_info, info_df, 
            by = "info_terms",
            all = TRUE
          )
        }

      }
      
      # add rownames and transpose additional_info:
      additional_info <- additional_info %>%
        column_to_rownames("info_terms")
      additional_info <- as.data.frame(t(additional_info))

      # convert to GRanges and remove alternate chromosome scaffolds:
      gr <- GRanges(
        seqnames = vcf_df$V1,
        ranges = IRanges(start = vcf_df$V2, vcf_df$V2),
        strand = "*",
        id = vcf_df$V3,
        join_base = vcf_df$V4,
        join = vcf_df$V5,
        join_chr = gsub(
          ":.*$", "", 
          gsub("^.*chr", "chr", vcf_df$V5)
        ),
        join_coord = as.numeric(
          gsub(
            "[^0-9.-]", "", 
            gsub("^.*chr.*:", "", vcf_df$V5)
          )
        ),
        quality = vcf_df$V6,
        filter = vcf_df$V7,
        DISC_MAPQ = additional_info$DISC_MAPQ,
        EVDNC = additional_info$EVDNC,
        type = additional_info$SVTYPE,
        MAPQ = additional_info$MAPQ,
        MATEID = additional_info$MATEID,
        MATENM = additional_info$MATENM ,
        NM = additional_info$NM,
        NUMPARTS = additional_info$NUMPARTS
      )

      return(gr)

    } else {
      print(paste0("VCF data is not in correct format for ", samplename,
        ", returning NA"))
      return(NA)
    }

  } else {
    print(paste0("VCF file empty or no VCF file exists for ", samplename,
        ", returning NA"))
    return(NA)
  }

  }