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
#' @examples
#' # Assuming 'se' is your sample
#' # Also assuming you have a valid output directory path in 'outputPath'
#' generateReport(se = se, tcgaData = tcgaData, outputDir = outputPath)
#'
#' @export
generateReport <- function(se, outputDir) {

  # Compute TCGA stats
  data("tcgaStats", package = "OvRSeq")

  # Iterate over each patient/sample in se
  patientIDs <- rownames(colData(se))
  for (patientID in patientIDs) {
    # ... (previous code to prepare patientData and reportData)

    # Convert reportData to a LaTeX table string
    latexTable <- print(xtable(reportData), include.rownames = FALSE, include.colnames = TRUE, only.contents = TRUE, comment = FALSE, sanitize.text.function = function(x) x)

    # Create a temporary Rmd file for each patient
    rmdFile <- tempfile(fileext = ".Rmd")

    # Write Rmd content with the LaTeX table directly embedded
    writeLines(con = rmdFile, text = c(
      "---",
      paste("title: 'OvRSeq Analysis Report for", patientID, "'"),
      "output: pdf_document",
      "---",
      "",
      "## OvRSeq Analysis Results for Patient:", patientID,
      "",
      "This report shows the analysis results for the patient using the OvRSeq pipeline, alongside TCGA reference statistics.",
      "",
      "### Patient Values and TCGA Reference",
      "",
      latexTable,  # Embed the LaTeX table directly
      ""
      # Additional content, plots, or tables can be added here
    ))

    # Render the Rmd file to PDF
    outputFilePath <- file.path(outputDir, paste0("OvRSeqReport_", patientID, ".pdf"))
    rmarkdown::render(input = rmdFile, output_file = outputFilePath)
  }
}
