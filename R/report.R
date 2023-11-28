#' Load TCGA Statistics Data
#'
#' This function loads the TCGA statistics data from the 'OvRSeq' package. The data is typically used for comparative analysis in the context of high-grade serous ovarian cancer (HGSOC) research.
#'
#' @return A data frame containing TCGA statistics.
#'
#' @details
#' The function retrieves the 'tcgaStats' dataset, which includes various statistical measures derived from The Cancer Genome Atlas (TCGA) data. This dataset is crucial for analyzing and comparing patient-specific data against a broader cancer database, facilitating insights into molecular patterns, biomarker distributions, and other critical aspects in HGSOC studies.
#'
#' @examples
#' tcga_stats <- load_tcgaStats()
#'
#' @export
load_tcgaStats <- function() {
  data("tcgaStats", package = "OvRSeq")
  return(tcgaStats)

}

#' Generate OvRSeq Analysis Reports for Patients
#'
#' This function generates individual PDF reports for each patient in a given dataset.
#' Each report includes the analysis results for the patient using the OvRSeq pipeline,
#' alongside TCGA (The Cancer Genome Atlas) reference statistics.
#'
#' @param se An object containing the patient sample data.
#' @param outputDir A string specifying the directory where the PDF reports will be saved.
#'
#' @return This function does not return a value, but generates PDF reports and saves them
#' in the specified output directory.
#'
#'@importFrom knitr kable
#'@importFrom ggplot2 ggsave
#'@importFrom stringr str_remove
#'@import dplyr
#'
#' @examples
#' # Assuming 'se' is your sample
#' # Also assuming you have a valid output directory path in 'outputPath'
#' generateReport(se = se, tcgaData = tcgaData, outputDir = outputPath)
#'
#' @export
OvRSeqReport <- function(se, outputDir) {
  # Compute TCGA stats
  tcgaStats <- load_tcgaStats()
  tcgaStats$IQR <- sprintf("%.2f (%.2f-%.2f)",
                           tcgaStats$`Median.50%`,
                           tcgaStats$`Q1.25%`,
                           tcgaStats$`Q3.75%`)

  # Ensure the output directory exists
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }

  escapeUnderscores <- function(text) {
    gsub("_", " ", text, fixed = TRUE)
  }
  escapePipe <- function(text) {
    gsub("|", "\\\\| ", text, fixed = TRUE)
  }
  # Iterate over each patient/sample in se
  patientIDs <- rownames(colData(se))
  for (patientID in patientIDs) {
    cat(paste0("[+] Creating report for: ", patientID))
    # Extract patient data
    patientColData <- colData(se)[patientID,]
    vector_genes <- c("CD274", "GZMB", "PRF1", "C1QA","CD8A", "IDO1", "FOXP3",
                      "TREM2", "STAT1", "HLA-DRA", "CXCL10")
    patientExpressionData <- t(assay(se)[vector_genes,patientID])
    # Keep only numeric columns
    patientData <- cbind(patientColData,patientExpressionData)
    #PLot
    p1 <- plot_vulnerabilitymap(se[,patientID])
    p2 <- plot_ggmarginal_sample(se[,patientID],color_var = "BRCAness")
    p3 <- plot_immune_signature_one_sample(se, patientID)
    p4 <- plot_quantiseq_one_sample(se, patientID)
    # Create a temporary file for the plot
    plotFile1 <- tempfile(fileext = ".jpeg")
    plotFile2 <- tempfile(fileext = ".jpeg")
    plotFile3 <- tempfile(fileext = ".jpeg")
    plotFile4 <- tempfile(fileext = ".jpeg")

    ggsave(plotFile1, plot = p1, width = 148, height = 100, units = "mm", dpi = 600)
    ggsave(plotFile2, plot = p2, width = 120, height = 100, units = "mm", dpi = 600)
    ggsave(plotFile3, plot = p3, width = 148, height = 100, units = "mm", dpi = 600)
    ggsave(plotFile4, plot = p4, width = 148, height = 100, units = "mm", dpi = 600)

    # Prepare patient data as a LaTeX table
    patientTable <- paste(
      "\\begin{tabular}{ll}",  # Define the table with 2 columns
      "\\hline",  # Horizontal line
      "Feature & Value \\\\",  # Table header
      "\\hline",  # Horizontal line
      sprintf("BRCAness status & %s (%.2f) \\\\", patientData$BRCAness, patientData$BRCAness_Prob),
      #"\\hline",  # Horizontal line
      paste0("Infiltration status & ", patientData$InfiltrationStatus,  " \\\\"),  # Row for Infiltration Status
      #"\\hline",  # Horizontal line
      paste0("Molecular subtypes & ", str_remove(escapeUnderscores(patientData$Tumor_Molecular_Subtypes), "consensus"),  " \\\\"),  # Row for Tumor Molecular Subtypes
      paste0("BRCAness immunotype & ", patientData$BRCAness_immunotype,  " \\\\"),  # Row for Tumor Molecular Subtypes
      sprintf("Vulnerability score & %.2f \\\\", patientData$Vulnerability_Score),
      paste0("Immuno phenoscore & ", patientData$IPS,  " \\\\"),  # Row for Tumor Molecular Subtypes
      sprintf("CYT to C1QA ratio (C2C) & %.2f \\\\", patientData$mapped_ratio_CYT_C1QA),
      paste0("Angiogenesis Score & ", "",  " \\\\"),  # Row for Tumor Molecular Subtypes

      "\\hline",  # Horizontal line
      "\\end{tabular}",
      sep = "\n")

    referenceTable <- paste(
      "\\begin{tabular}{lll}",  # Define the table with 3 columns
      "\\hline",  # Horizontal line
      "Feature & Patient Value & TCGA IQR \\\\",  # Table header
      "\\hline",  # Horizontal line
      paste(sapply(vector_genes, function(gene) {
        sprintf("%s & %.2f & %s \\\\", gene, patientData[[gene]], tcgaStats[gene, "IQR"])
      }
      ), collapse = "\n"),
      "\\hline",  # Horizontal line
      "\\end{tabular}",
      sep = "\n"
    )


    # Create a temporary Rmd file for each patient
    rmdFile <- tempfile(fileext = ".Rmd")
    # Write Rmd content
    writeLines(con = rmdFile, text = c(
      "---",
      "output: pdf_document",
      "geometry: margin = 0.6in",
      "header-includes:",
      "   - \\usepackage{multicol}",
      "   - \\usepackage{graphicx}",
      "---",
      "",
      paste0("## OvRSeq Analysis Report for ", patientID),
      "",
      as.character(Sys.Date()),
      "\\vspace{5mm}",
      "",
      "The vulnerability map indicate based on BRCAness probability and CYT to C1QA ratio (C2C) indications with a high vulnerability (score) for response to combination immunotherapy with PARPi and immune checkpoint inhibitors.",
      "",
      "\\begin{multicols}{2}",  # Start first two-column layout
      "",
      paste0("\\includegraphics{", gsub("\n", "", plotFile1), "}"),  # Include the first plot with newline removed
      "",
      "\\columnbreak",  # Break to the second column
      "",
      "\\vspace{10mm}",
      "\\textbf{Patient values}",
      "\\vspace{5mm}",
      "",
      patientTable,  # Add patient information here
      "",
      "\\end{multicols}",  # End first two-column layout
      "",
      "\\textbf{Molecular markers and TCGA-OV Reference Values}. Marker gene expression levels and respective Q1 and Q3 levels (interquartile range) from the TCGA cohort.",
      "",
      "\\begin{multicols}{2}",  # Start second two-column layout
      "",
      paste0("\\includegraphics{", gsub("\n", "", plotFile2), "}"),  # Include the second plot with newline removed
      "",
      "\\columnbreak",  # Break to the second column in the second layout
      "",
      referenceTable,  # Add reference table here
      "",
      "\\end{multicols}",  # End second two-column layout
      "",
      "\\textbf{Marker immune signatures and pathway scores and quantified immune cell infiltrates}, integrating GSVA-derived metrics and quanTIseq estimations for a detailed immunological profile.",
      "",
      "\\begin{multicols}{2}",  # Start first two-column layout
      "",
      paste0("\\includegraphics{", gsub("\n", "", plotFile3), "}"),  # Include the first plot with newline removed
      "",
      "\\columnbreak",  # Break to the second column
      "",
      paste0("\\includegraphics{", gsub("\n", "", plotFile4), "}"),  # Include the first plot with newline removed
      "",
      "\\end{multicols}"  # End first two-column layout
    ))

    # Render the Rmd file to PDF
    outputFilePath <- file.path(outputDir, paste0("OvRSeqReport_", patientID, ".pdf"))
    rmarkdown::render(input = rmdFile, output_file = outputFilePath)

    # Optionally, delete the temporary plot file
    unlink(plotFile1)
    unlink(plotFile2)
    unlink(plotFile3)
    unlink(plotFile4)

  }
}
