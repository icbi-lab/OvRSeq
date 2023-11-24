calculateTcgaStats <- function(tcgaData) {
  if (!inherits(tcgaData, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }

  # Extract colData
  TCGA_OV <- OvRSeq(TCGA_OV, normalize = F)
  tcgaColData <- colData(TCGA_OV)
  tcgaExpressionData <- t(assay(TCGA_OV)[c("CD274", "GZMB", "PRF1", "C1QA","CD8A"),])
  # Keep only numeric columns
  tcgaColData <- cbind(tcgaColData,tcgaExpressionData)
  numericCols <- sapply(tcgaColData, is.numeric)
  tcgaColData <- tcgaColData[, numericCols]

  # Calculate median, 25th percentile (Q1), and 75th percentile (Q3) for each column
  stats <- t(apply(tcgaColData, 2, function(x) {
    quantiles <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    c(Q1 = round(quantiles[1],3), Median = round(quantiles[2],3), Q3 = round(quantiles[3],3))
  }))

  # Convert to data frame
  statsDf <- as.data.frame(stats)
  rownames(statsDf) <- NULL
  statsDf$Feature <- colnames(tcgaColData)
  rownames(statsDf) <- statsDf$Feature
  tcgaStats <- statsDf
  usethis::use_data(tcgaStats, overwrite = T)
  return(statsDf)
}
