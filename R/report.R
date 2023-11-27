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
#'@import dplyr
#'
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
  tcgaStats$IQR <- paste0(tcgaStats$`Median.50%`," (",tcgaStats$`Q1.25%`,"-",tcgaStats$`Q3.75%`,")")
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
    # Create a temporary file for the plot
    plotFile1 <- tempfile(fileext = ".jpeg")
    plotFile2 <- tempfile(fileext = ".jpeg")

    ggsave(plotFile1, plot = p1, width = 148, height = 148, units = "mm", dpi = 600)
    ggsave(plotFile2, plot = p2, width = 148, height = 148, units = "mm", dpi = 600)

    # Prepare patient data as a LaTeX table
    patientTable <- paste(
      "\\begin{tabular}{cc}",  # Define the table with 2 columns
      "\\hline",  # Horizontal line
      "Feature & Value \\\\",  # Table header
      "\\hline",  # Horizontal line
      paste0("BRCAness Status & ", patientData$BRCAness," (", patientData$BRCAness_Prob,")", " \\\\"),  # Row for BRCAness Status
      #"\\hline",  # Horizontal line
      paste0("Infiltration Status & ", patientData$InfiltrationStatus,  " \\\\"),  # Row for Infiltration Status
      #"\\hline",  # Horizontal line
      paste0("Molecular Subtypes & ", escapeUnderscores(patientData$Tumor_Molecular_Subtypes),  " \\\\"),  # Row for Tumor Molecular Subtypes
      paste0("BRCAness immunotype & ", patientData$BRCAness_immunotype,  " \\\\"),  # Row for Tumor Molecular Subtypes
      paste0("Vulnerability Score & ", round(patientData$Vulnerability_Score,3),  " \\\\"),  # Row for Tumor Molecular Subtypes
      paste0("Immuno Phenoscore & ", round(patientData$IPS,3),  " \\\\"),  # Row for Tumor Molecular Subtypes
      paste0("CYT to C1QA Ratio & ", round(patientData$mapped_ratio_CYT_C1QA,3),  " \\\\"),  # Row for Tumor Molecular Subtypes
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
          paste0(gene, " & ", round(patientData[[gene]], 3), " & ", tcgaStats[gene, "IQR"], " \\\\")
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
      "title: 'OvRSeq Analysis Report for ", patientID, "'",
      "output: pdf_document",
      'date: "`r Sys.Date()`"',
      "header-includes:",
      "   - \\usepackage{multicol}",
      "   - \\usepackage{graphicx}",
      "---",
      "",
      paste0("## OvRSeq Analysis Results for Patient: ", patientID),
      "",
      "The report starts with a Vulnerability Map highlighting HGSOC aspects, including BRCAness Status, Infiltration Status, and Tumor Molecular Subtypes, setting the stage for a personalized molecular analysis.",
      "",
      "\\begin{multicols}{2}",  # Start first two-column layout
      "",
      paste0("\\includegraphics{", gsub("\n", "", plotFile1), "}"),  # Include the first plot with newline removed
      "",
      "\\columnbreak",  # Break to the second column
      "",
      "\\textbf{Patient Values}",
      "",
      patientTable,  # Add patient information here
      "",
      "\\end{multicols}",  # End first two-column layout
      "",
      "\\textbf{Molecular markers and TCGA-OV Reference Values}",
      "",
      "This report section compares the patient's molecular markers with TCGA interquartile ranges (IQR). The table includes key markers like Immuno Phenoscore and CYT to C1QA ratio, crucial for personalized cancer profiling and treatment.",
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
      ""
      # ... rest of the report content
    ))



    # Render the Rmd file to PDF
    outputFilePath <- file.path(outputDir, paste0("OvRSeqReport_", patientID, ".pdf"))
    rmarkdown::render(input = rmdFile, output_file = outputFilePath)

    # Optionally, delete the temporary plot file
    unlink(plotFile1)
    unlink(plotFile2)
  }
}

# ... following code ...

