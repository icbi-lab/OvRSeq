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
  data("tcgaStats", package = "OvRSeq")
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
    patientExpressionData <- t(assay(se)[c("CD274", "GZMB", "PRF1", "C1QA","CD8A"),patientID])
    # Keep only numeric columns
    patientData <- cbind(patientColData,patientExpressionData)
    #PLot
    p1 <- plot_vulnerabilitymap(se[,patientID])
    # Create a temporary file for the plot
    plotFile1 <- tempfile(fileext = ".jpeg")

    ggsave(plotFile1, plot = p1, width = 148, height = 148, units = "mm", dpi = 600)
    # Extract common features and filter patient data

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
      paste0("Tumor Molecular Subtypes & ", escapeUnderscores(patientData$Tumor_Molecular_Subtypes),  " \\\\"),  # Row for Tumor Molecular Subtypes
      paste0("BRIT & ", "",  " \\\\"),  # Row for Tumor Molecular Subtypes
      paste0("Angiogenesis Score & ", "",  " \\\\"),  # Row for Tumor Molecular Subtypes
      "\\hline",  # Horizontal line
      "\\end{tabular}",
      sep = "\n")

    referenceTable <- paste(
      "\\begin{table}[h!]",
      "\\centering",
      "\\begin{tabular}{lll}",  # Define the table with 3 columns
      "\\hline",  # Horizontal line
      "Feature & Patient Value & TCGA IQR \\\\",  # Table header
      "\\hline",  # Horizontal line
      paste0("Immuno Phenoscore & ", patientData$IPS, " & ", tcgaStats["IPS", "IQR"], " \\\\"),  # Row for Immuno Phenoscore
      paste0("CYT to C1QA Ratio & ", patientData$ratio_CYT_C1QA, " & ", tcgaStats["ratio_CYT_C1QA", "IQR"], " \\\\"),  # Row for CYT to C1QA Ratio
      paste0("IFNG Higgs & ", round(patientData$`IFNG Higgs`,3), " & ", tcgaStats["IFNG Higgs", "IQR"], " \\\\"),  # Row for IFNG Higgs
      paste0("IFNG Ayers & ", round(patientData$`IFNG Ayers`,3), " & ", tcgaStats["IFNG Ayers", "IQR"], " \\\\"),  # Row for IFNG Ayers
      paste0("CD274 & ", round(patientData$CD274,3), " & ", tcgaStats["CD274", "IQR"], " \\\\"),  # Row for CD274
      paste0("GZMB & ", round(patientData$GZMB,3), " & ", tcgaStats["GZMB", "IQR"], " \\\\"),  # Row for GZMB
      paste0("PRF1 & ", round(patientData$PRF1,3), " & ", tcgaStats["PRF1", "IQR"], " \\\\"),  # Row for PRF1
      paste0("C1QA & ", round(patientData$C1QA,3), " & ", tcgaStats["C1QA", "IQR"], " \\\\"),  # Row for C1QA
      paste0("CD8A & ", round(patientData$CD8A,3), " & ", tcgaStats["CD8A", "IQR"], " \\\\"),  # Row for CD8A
      "\\hline",  # Horizontal line
      "\\end{tabular}",
      "\\end{table}",
      sep = "\n")


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
      "## OvRSeq Analysis Results for Patient: ", patientID,
      "",
      "The report begins with the Vulnerability Map, illustrating key aspects of high-grade serous ovarian cancer (HGSOC). It includes the BRCAness Status, showing the probability of BRCA-related vulnerabilities, Infiltration Status, indicating immune system engagement, and Tumor Molecular Subtypes, essential for targeted therapeutic strategies. This section lays the groundwork for a detailed patient-specific molecular analysis.",
      "",
      "\\begin{multicols}{2}",  # Start two-column layout
      "",
      paste0("\\includegraphics{", gsub("\n", "", plotFile1), "}"),  # Include the plot with newline removed
      "",
      "\\columnbreak",  # Break to the second column
      "",
      "\\textbf{Patient Values and TCGA Reference}",
      "",
      # Add patient information here
      patientTable,
      "",
      "\\end{multicols}",  # End two-column layout
      "",
      "\\textbf{Molecular markers and TCGA-OV Reference Values}",

      # ... rest of the report content
      "This report section compares the patient's molecular markers with TCGA interquartile ranges (IQR). The table includes key markers like Immuno Phenoscore and CYT to C1QA ratio, crucial for personalized cancer profiling and treatment.",
      referenceTable
    ))


    # Render the Rmd file to PDF
    outputFilePath <- file.path(outputDir, paste0("OvRSeqReport_", patientID, ".pdf"))
    rmarkdown::render(input = rmdFile, output_file = outputFilePath)

    # Optionally, delete the temporary plot file
    unlink(plotFile1)
  }
}
