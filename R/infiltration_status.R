#' Load Tumor Immune Phenotype Signature
#'
#' Loads the tumor immune phenotype signature from the
#' "tumor_immune_phenotype_signature.rda" data file. This signature contains
#' a list of genes associated with the immune phenotype of tumors.
#'
#' @return
#' A character vector of gene symbols representing the tumor immune phenotype signature.
#'
#' @usage
#' tumor_immune_phenotype_signature <- load_tumor_immune_phenotype_signature()
#'
#' @examples
#' # Load the tumor immune phenotype signature
#' tumor_immune_phenotype_signature <- load_tumor_immune_phenotype_signature()
#'
#' @export
load_tumor_immune_phenotype_signature <- function() {
  # Load the data file
  data("tumor_immune_phenotype_signature",package = "OvRSeq")
  return(tumor_immune_phenotype_signature)
}



#' Train Random Forest Classifier for Tumor Infiltration Status
#'
#' Trains a random forest classifier using a predefined gene signature to predict tumor immune phenotypes,
#' such as infiltrate, excluded, or desert based on the gene expression profiles from a SummarizedExperiment object.
#'
#' @param se A SummarizedExperiment object containing gene expression data and a column `TinfStatus` in its `colData`
#' which indicates the tumor immune phenotype.
#' @param gene_signature A character vector containing the names of genes in the signature to be used for model training.
#' @param num_estimators The number of trees to grow in the random forest model. Default is 300.
#'
#' @return A random forest model object trained on the gene expression data and tumor immune phenotype status.
#'
#' @examples
#' # Assuming 'se' is a preloaded SummarizedExperiment object with the required data
#' # and 'tumor_immune_phenotype_signature' is the predefined gene signature list:
#' rf_classifier <- train_rf_classifier_infiltration_status(se, tumor_immune_phenotype_signature)
#'
#' @details
#' The function extracts the expression data for the genes in the provided signature from a SummarizedExperiment object
#' and uses the `TinfStatus` from the column metadata to train a random forest classifier.
#' The `TinfStatus` should be a factor that represents the tumor immune phenotype status with levels such as
#' 'infiltrate', 'excluded', or 'desert'. The gene expression matrix is transposed to ensure proper dimensions
#' for the random forest training function. The `train` function from the `caret` package is used for training
#' the model, with the number of trees set by the `num_estimators` parameter.
#'
#' The function sets the seed for reproducibility and ensures the output is a trained model
#' that can predict the `TinfStatus` based on gene expression profiles.
#'
#' @importFrom caret train
#' @references
#' Desbois M, Udyavar AR, Ryner L, Kozlowski C, Guan Y, Dürrbaum M, et al. Integrated digital
#' pathology and transcriptome analysis identifies molecular mediators of T-cell exclusion in
#' ovarian cancer. Nat Commun. 2020;11:5583.
#'
#' @export
train_rf_classifier_infiltration_status <- function(se, gene_signature, num_estimators = 300) {
  set.seed(42) # Ensure reproducibility

  # Extract the gene expression data and the response variable 'TinfStatus'
  X <- t(assay(se))
  y <- colData(se)$TinfStatus

  # Ensure the gene signature is present in the data
  if (!all(gene_signature %in% rownames(se))) {
    stop("Not all genes in the signature are present in the expression dataset.")
  }

  # Subset the gene expression matrix to the signature genes
  X <- X[, gene_signature]

  # Configure the control parameters for the training process
  train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)

  # Train the Random Forest Classifier
  rf_model <- train(X, as.factor(y), method = "rf", trControl = train_control,
                    metric = "Accuracy", ntree = num_estimators)

  # Return the trained classifier
  return(rf_model)
}



#' Load Infiltration Status Classifier
#'
#' Loads the pre-trained classifier for tumor immune infiltration status from the
#' "classifier_infiltration_status.rda" data file. This classifier is used to determine
#' the infiltration status (e.g., infiltrated, excluded, desert) within the tumor microenvironment.
#'
#' @return
#' A pre-trained classifier object. Typically, this is a machine learning model object
#' that can be used to predict infiltration status based on gene expression profiles.
#'
#' @usage
#' classifier_infiltration_status <- load_classifier_infiltration_status()
#'
#' @examples
#' # Load the Infiltration Status classifier
#' classifier_infiltration_status <- load_classifier_infiltration_status()
#'
#' # Now you can use `classifier_infiltration_status` with gene expression data to
#' # predict infiltration statuses.
#'
#' @export
load_classifier_infiltration_status <- function() {
  data("classifier_infiltration_status",package = "OvRSeq")
  return(classifier_infiltration_status)
}


#' Classify Tumor Immune Infiltration Status
#'
#' This function applies a trained classifier to a SummarizedExperiment object to classify
#' samples based on their tumor immune infiltration status using a specified gene signature.
#'
#' @param se A SummarizedExperiment object containing gene expression data where samples are
#'   columns and genes are rows.
#' @param classifier A trained classifier object capable of making predictions based on the
#'   gene expression data provided in `se`.
#' @param gene_signature A character vector containing gene symbols that make up the gene
#'   signature used by the classifier. The function will subset the expression data in `se`
#'   based on these genes to make predictions.
#'
#' @return The SummarizedExperiment object `se` with an additional column in `colData` named
#'   `InfiltrationStatus`, which contains the classification results for each sample.
#'
#' @details
#' The function first validates that the provided `se` object is indeed a SummarizedExperiment.
#' It then extracts the expression data for the genes in the `gene_signature` from the `se` object.
#' This gene expression matrix is transposed (samples as rows, genes as columns) and used as input
#' for the classifier to predict the infiltration status. The predictions are then added to the
#' `colData` of the `se` object, allowing for easy retrieval and further analysis.
#'
#' @examples
#' # Assuming `se` is a SummarizedExperiment object with gene expression data,
#' # `rf_classifier` is a trained random forest classifier, and
#' # `immune_genes` is a vector of gene symbols making up the immune signature:
#' se <- classify_infiltration_status(se, rf_classifier, immune_genes)
#'
#' @export
classify_infiltration_status <- function(se, classifier, gene_signature) {
  # Check if 'se' is a SummarizedExperiment object
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }

  # Check gene_signature
  gene_signature <- intersect(rownames(se), gene_signature)
  # Extract gene expression matrix for the gene signature
  gene_expression <- assay(se)[gene_signature, , drop = FALSE]

  # Predict infiltration status using the classifier
  predicted_status <- predict(classifier, t(gene_expression))

  # Add the predictions to colData
  colData(se)$InfiltrationStatus <- predicted_status

  # Print message to console
  cat("Classification complete: ", length(pred), "samples were classified for Infiltration status.\n")
  # Return the updated SummarizedExperiment object
  return(se)
}

