#' Train a random forest model on gene expression data
#'
#' This function trains a random forest model on a summarized experiment object
#' containing gene expression data, using a specified gene signature to select a
#' subset of genes, and a specified label in the colData to use as the target variable.
#' The function returns the best random forest model.
#'
#' @param se A SummarizedExperiment object containing gene expression data
#' @param label A string indicating the name of the column in colData containing the label to use as the target variable for classification
#' @param gene_signature A character vector containing the names of genes to use in the analysis
#' @return A trained random forest model
#' @importFrom caret train
#' @examples
#' # Load the TCGA_OV.rda and BRCAness_signature.rda files
#' data(TCGA_OV)
#' data("brcaness_signature")
#'
#' # Train a random forest model using the BRCAness label
#' rf_model <- train_rf(se, "BRCAness", brcaness_signature)
#'
#' # Print the model
#' print(rf_model)
#'
#' # Make predictions using the model
#' predictions <- predict(rf_model, newdata = t(assay(TCGA_OV)[rownames(BRCAness_signature),]))
#'
#' # Print the predictions
#' print(predictions)
#' @export
train_rf <- function(se, label,gene_signature){

  # Subset data to intersecting genes
  gene_names <- intersect(rownames(se), gene_signature)
  count_data_subset <- assay(se)[gene_names, ]

  # Subset BRCAness label and remove missing values
  sample_label <- colData(se)[,label]
  sample_label <- sample_label[!is.na(sample_label)]
  sample_label <- factor(unlist(lapply(sample_label, function(x){ifelse(x==sample_label[1],1,0)})), levels = c(0,1))

  count_data_subset <- count_data_subset[, !is.na(sample_label)]

  set.seed(42)
  rf_model <- train(t(count_data_subset), sample_label,method = "rf", metric= "Accuracy")

  return(rf_model)
}




#' Load BRCAness Classifier
#'
#' Loads the pre-trained BRCAness classifier from the "brcaness_classifier.rda"
#' data file and returns it as the `brcaness_classifier` object.
#'
#' @return
#' The pre-trained BRCAness classifier as a list.
#'
#' @usage
#' load_brcaness_classifier()
#'
#' @examples
#' # Load the BRCAness classifier
#' brcaness_classifier <- load_brcaness_classifier()
#'
#' @export
load_brcaness_classifier <- function() {
  data("brcaness_classifier",package = "OvRSeq")
  return(brcaness_classifier)

}



#' Classify Samples Using BRCAness Classifier
#'
#' This function uses the BRCAness classifier to predict the BRCAness status of
#' the samples in a SummarizedExperiment object. It first checks that all the
#' genes in the BRCAness signature are present in the count data of the
#' SummarizedExperiment object before making predictions.
#'
#' @param se A SummarizedExperiment object with count data.
#' @param brcaness_classifier A random forest classifier object trained on
#'   BRCAness samples.
#' @param brcaness_signature A vector of gene symbols representing the BRCAness
#'   gene signature.
#' @return A SummarizedExperiment with BRCAness predictions, in ColData.
#' @export
classify_brcaness <- function(se, brcaness_classifier, brcaness_signature) {

  # Check that 'se' is a SummarizedExperiment object
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }

  # Check that all genes in BRCAness signature are present in count data
  if (!all(brcaness_signature %in% rownames(assay(se)))) {
    stop("Not all genes in BRCAness signature are present in the count data.")
  }

  # Subset count data to BRCAness signature genes
  count_data_subset <- assay(se)[brcaness_signature, ]

  # Check for the presence of genes in the subset
  if (nrow(count_data_subset) != length(brcaness_signature)) {
    stop("Gene signature subset did not return the expected number of genes.")
  }

  # Make predictions
  pred <- predict(brcaness_classifier, t(count_data_subset))
  colData(se)$BRCAness <- pred


  # Get prediction probabilities
  # The exact method may vary depending on the classifier used
  predProbs <- predict(brcaness_classifier, t(count_data_subset), type = "prob")[,"1"]
  colData(se)$BRCAness_Prob <- predProbs

  # Print message to console
  cat("Classification complete: ", length(pred), "samples were classified for BRCAness status.\n")

  return(se)
}


#' Stratify HGSOC Patients Based on BRCAness and Immune Subtypes
#'
#' This function stratifies patients in a `SummarizedExperiment` object into categories based on BRCAness status, tumor immunephenotype, and molecular subtype. It identifies patients with a BRCAness immunetype (BRIT), which combines the BRCAness phenotype with an infiltrated tumor immune phenotype and an immune-reactive molecular subtype (IMR).
#'
#' @param se A `SummarizedExperiment` object containing patient data, including BRCAness status, tumor immunephenotype, and molecular subtype.
#'
#' @return The modified `SummarizedExperiment` object with an added `immunotype` column in `colData`, classifying each patient as 'BRIT', 'noBRIT', or 'other'.
#'
#' @details
#' The function uses a list of 157 genes and random forest analysis to classify the tumor immunephenotype of each patient as infiltrated, excluded, or desert. The BRIT group consists of patients with BRCAness, an infiltrated tumor immune phenotype, and an immune-reactive molecular subtype (IMR). This group is potentially more responsive to PARP inhibitor and immune checkpoint inhibitor therapies. The function also accounts for the presence of suppressive immune cells and does not show a significant difference in overall survival between BRIT and noBRIT groups. This stratification provides insights into the complex interplay between genetic predisposition and immune response in HGSOC, aiding in personalized treatment planning.
#'
#' @examples
#' # se is a pre-loaded SummarizedExperiment object
#' se <- BRCAness_immunotype(se)
#'
#' @importFrom SummarizedExperiment colData
#' @export
BRCAness_immunotype <- function(se){
  # Check and compute necessary values if they are not present
  if (!("BRCAness" %in% colnames(colData(se)))) {
    # Compute BRCAness here (placeholder, replace with actual computation)
    # BRCAness classification
    tryCatch({
      brcaness_classifier <- load_brcaness_classifier()
      brcaness_signature <- load_brcaness_signature()
      se <- classify_brcaness(se, brcaness_classifier, brcaness_signature)
    }, warning = function(w) {
      warning("Failed to classify BRCAness: ", w$message)
    }, error = function(e) {
      warning("Error in classifying BRCAness: ", e$message)
    })
  }
  if (!("InfiltrationStatus" %in% colnames(colData(se)))) {
    ## Infiltration status classification
    tryCatch({
      classifier_infiltration_status <- load_classifier_infiltration_status()
      immune_phenotype_signature <- load_tumor_immune_phenotype_signature()
      se <- classify_infiltration_status(se, classifier_infiltration_status, immune_phenotype_signature)
    }, warning = function(w) {
      warning("Failed to classify infiltration status: ", w$message)
    }, error = function(e) {
      warning("Error in classifying infiltration status: ", e$message)
    })
  }
  if (!("Tumor_Molecular_Subtypes" %in% colnames(colData(se)))) {
    # Consensus molecular subtypes
    tryCatch({
      se <- get_consensus_ov_subtypes(se, ids_type = "symbol")
    }, warning = function(w) {
      warning("Failed to determine consensus molecular subtypes: ", w$message)
    }, error = function(e) {
      warning("Error in determining consensus molecular subtypes: ", e$message)
    })
  }

  # Compute immunotype based on existing or computed data
  data <- colData(se)
  immunotype <- apply(data, 1, function(sample){
    if (sample["BRCAness"] == 1 && sample["InfiltrationStatus"] == "Infiltrated" && sample["Tumor_Molecular_Subtypes"] == "IMR_consensus") {
      return("BRIT")
    } else if (sample["BRCAness"] == 1 && sample["InfiltrationStatus"] != "Infiltrated" && sample["Tumor_Molecular_Subtypes"] != "IMR_consensus") {
      return("noBRIT")
    } else {
      return("other")
    }
  })
  colData(se)$BRCAness_immunotype <- immunotype
  cat("Classification complete: ", length(pred), "samples were classified for BRCAness immunotype\n")
  return(se)
}


