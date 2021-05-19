load_breakpoints <- function(samplenames, in_path, high_conf = TRUE) {

    # load VCFs and convert to granges:
    breakpoint_gr <- lapply(samplenames, function(x) {
    
      print(paste0("Loading ", x, "..."))
    
      if (high_conf) {
        in_file <- paste0(in_path, x, "/", x, ".svaba.sv.vcf")
      } else {
        in_file <- paste0(in_path, x, "/", x, ".svaba.semifiltered.sv.formatted.vcf")
      }

      vcf_df <- tryCatch(
        read.table(
          in_file,
          sep = "\t"
        ),
        warning = function(w) {NA},
        error = function(e) {NA}
      )
      
      if (!is.na(vcf_df)) {
  
        # remove alternate chromosomes:
        vcf_df <- vcf_df[grep("chrG|chrK|chrJ", vcf_df$V1, invert = T),]
        vcf_df <- vcf_df[grep("chrG|chrK|chrJ", vcf_df$V5, invert = T),]

        if (is.data.frame(vcf_df)) {

          # split additional info field:
          additional_info <- gsub("^.*=", "", do.call("rbind", strsplit(vcf_df$V8, ";")))
    
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
            DISC_MAPQ = additional_info[,1],
            EVDNC = additional_info[,2],
            type = additional_info[,3],
            MAPQ = additional_info[,4],
            MATEID = additional_info[,5],
            MATENM = additional_info[,6],
            NM = additional_info[,7],
            NUMPARTS = additional_info[,8]
          )
    
          return(gr)

        } else {
          return(NA)
        }
    
        
    
      } else {
        return(NA)
      }
    
    })
    names(breakpoint_gr) <- samplenames

    return(breakpoint_gr)

  }