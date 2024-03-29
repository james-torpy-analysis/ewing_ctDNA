compare_VAF <- function(VAF_df, lab1, lab2, lim = 40, cortype = "pearson") {
  
  # remove samples without two values:
  VAF_df <- VAF_df[!is.na(VAF_df$VAF1) & !is.na(VAF_df$VAF2),]
  VAF_df <- VAF_df[VAF_df$VAF1 != "not_detected" & VAF_df$VAF2 != "not_detected",]
  VAF_df$VAF1 <- as.numeric(VAF_df$VAF1)
  VAF_df$VAF2 <- as.numeric(VAF_df$VAF2)
  
  # remove samples above lim:
  VAF_df <- VAF_df[VAF_df$VAF1 <= lim & VAF_df$VAF2 <= lim,]
  
  if (all(VAF_df$VAF1 < 1)) {
    VAF_df$VAF1 <- VAF_df$VAF1*100
  }
  if (all(VAF_df$VAF2 < 1)) {
    VAF_df$VAF2 <- VAF_df$VAF2*100
  }
  
  # calculate correlation between ddPCR and Andre VAFs:
  corr <- try(
    cor.test(x = VAF_df$VAF1, y = VAF_df$VAF2, method = cortype ), silent = T )
  
  # Fit regression line
  require(stats)
  reg <- try(lm(VAF2 ~ VAF1, data = VAF_df), silent = T)
  if (class(reg) != "try-error") {
    coeff=coefficients(reg)
  }

  p <- ggplot(
    VAF_df, aes(x = VAF1, y = VAF2, color = treatment) )
  p <- p + geom_point()
  p <- p + theme_cowplot(12)
  p <- p + xlim(c(0, lim))
  p <- p + ylim(c(0, lim))
  p <- p + xlab(paste0(lab1, " VAF"))
  p <- p + ylab(paste0(lab2, " VAF"))
  p <- p + geom_text_repel(data=VAF_df, aes(label=id), size = 3)
  if (class(corr) != "try-error") {
    p <- p + annotate("text", x = 20, y = 30, label = paste0(
      "R2=", round(corr$estimate, 2), ", p=", formatC(corr$p.value, format = "e", digits = 2) ), 
      color='red', size = 3.5 )
  }
  if (exists("coeff")) {
    p <- p + geom_abline(intercept = coeff[1], slope = coeff[2], color = "red")
  }
  return(p)
}