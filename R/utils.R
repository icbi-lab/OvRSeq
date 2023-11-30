calculateTcgaStats <- function(tcgaData) {
  if (!inherits(tcgaData, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }

  # Extract colData
  load_TCGA_OV()
  TCGA_OV <- OvRSeq(TCGA_OV, normalize = F)
  TCGA_OV <- computeCytC1qaRatio(TCGA_OV)
  tcgaColData <- colData(TCGA_OV)
  vector_genes <- c("CD274", "GZMB", "PRF1", "C1QA","CD8A", "IDO1", "FOXP3",
    "TREM2", "STAT1", "HLA-DRA", "CXCL10")
  tcgaExpressionData <- t(assay(TCGA_OV)[vector_genes,])
  # Keep only numeric columns
  tcgaColData <- cbind(tcgaColData,tcgaExpressionData)
  numericCols <- sapply(tcgaColData, is.numeric)
  tcgaColData <- tcgaColData[, numericCols]

  # Calculate median, 25th percentile (Q1), 75th percentile (Q3), minimum (Min), and maximum (Max) for each column
  stats <- t(apply(tcgaColData, 2, function(x) {
    quantiles <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    minVal <- min(x, na.rm = TRUE)
    maxVal <- max(x, na.rm = TRUE)
    c(Q1 = round(quantiles[1],3), Median = round(quantiles[2],3), Q3 = round(quantiles[3],3), Min = round(minVal,3), Max = round(maxVal,3))
  }))


  # Convert to data frame
  Feature <- rownames(stats)
  statsDf <- as.data.frame(stats)
  statsDf$Feature <-  Feature
  # Keep only unique rows based on the Feature column
  statsDf <- statsDf[!duplicated(statsDf$Feature),]
  rownames(statsDf) <- statsDf$Feature
  tcgaStats <- statsDf
  usethis::use_data(tcgaStats, overwrite = T)
  return(statsDf)
}



load_and_compute_thresholds <- function() {

  # Load the TCGA data
  se <- TCGA_OV
  data <- t(assay(se))

  # Extract survival columns
  SURV_TIME <- se$OS_MONTHS  # or the correct column name for survival time
  SURV_EVENT <- se$OS        # or the correct column name for survival event

  # Define the genes of interest
  genes_of_interest <- c("VEGFA", "KDR", "PDGFA", "PDGFB", "PDGFRA", "PDGFRB", "KIT")

  # Initialize a list to store thresholds
  thresholds <- list()

  # Compute the threshold for each gene using cutpointr
  for (gene in genes_of_interest) {
    if (gene %in% colnames(data)) {
      # Select the gene expression data and the survival data, removing missing values


      # Use cutpointr to compute the optimal threshold
      optimal_cut <- cutpointr(
        data[,gene],
        SURV_EVENT,
        pos_class = 1,
        neg_class = 0,
        method = maximize_metric,
        metric = youden,
        na.rm = T)  # adjust by parameter as needed


      # Store the threshold
      thresholds[[gene]] <- optimal_cut$optimal_cutpoint
    } else {
      warning(paste("Gene", gene, "not found in TCGA data."))
    }
  }

  # Return the list of thresholds
  return(thresholds)
}

