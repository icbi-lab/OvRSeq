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

train_rf_classifier_infiltration_status <- function(se, gene_signature, num_estimators = 300) {
  set.seed(42)
  # Extract gene expression matrix and TinfStatus from the SummarizedExperiment
  X <- t(assay(se))
  y <- colData(se)$TinfStatus

  # Subset the gene expression matrix based on the gene signature
  X <- X[, gene_signature]

  # Train the Random Forest Classifier
  rf_model <- train(X, as.factor(y),method = "rf", metric= "Accuracy")

  # Return the trained classifier
  return(rf_classifier)
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
  load("brcaness_classifier.rda")
  if (exists("brcaness_classifier", inherits = FALSE)) {
    return(brcaness_classifier)
  } else {
    stop("brcaness_classifier not found in the loaded data.")
  }
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

  # Check that all genes in BRCAness signature are present in count data
  if (!all(brcaness_signature %in% rownames(assay(se)))) {
    stop("Not all genes in BRCAness signature are present in the count data.")
  }

  # Subset count data to BRCAness signature genes
  count_data_subset <- assay(se)[brcaness_signature, ]

  # Make predictions
  pred <- predict(brcaness_classifier, t(count_data_subset))
  colData(se)$BRCAness <- pred

  return(se)
}

#' Obtain Consensus Subtypes of High-Grade Serous Ovarian Cancer
#'
#' This function gets the consensus molecular subtypes of high-grade serous (HGS) ovarian cancer using the `consensusOV` package.
#' The function implements a consensus classifier which consolidates and improves on the robustness of proposed subtype classifiers.
#'
#' @param se A Summarized Experiments of gene expression values, with gene symbol or entrez id as identifiers.
#' @param ids_type A character string indicating the type of IDs used in the 'data' row names. Can be either "symbol" or "entrez".
#' @param concordant.tumors.only Logical. If TRUE, only tumors that are concordantly classified across all datasets are included in the final classification. Defaults to TRUE.
#' @param remove.using.cutoff Logical. If TRUE, tumors with poor classification confidence are removed using the optimal cutoff. Defaults to FALSE.
#' @param percentage.dataset.removed Numeric value indicating the percentage of the dataset to be removed based on classification confidence. Defaults to 0.75.
#' @return Summarized Experiments with Tumor_Molecular_Subtypes of consensus molecular subtypes for the samples, in colData.
#' @importFrom consensusOV get.consensus.subtypes
#' @examples
#' load_TCGA_OV
#' TCGA_OV <- get_consensus_ov_subtypes(TCGA_OV, ids_type = "symbol")
#' @export
#' @references
#' Chen G, Kannan L, Geistlinger L, Kofia V, Safikhani Z, Gendoo D, Parmigiani G, Birrer M, Haibe-Kains B, Waldron L (2018).
#' "Consensus on molecular subtypes of high-grade serous ovarian carcinoma."
#' Clinical Cancer Research, 24, 4990. doi:10.1158/1078-0432.CCR-18-0784.
get_consensus_ov_subtypes <- function(se, ids_type,
                                      concordant.tumors.only = TRUE, remove.using.cutoff = FALSE,
                                      percentage.dataset.removed = 0.75) {

  if (ids_type == "symbol") {
    # Convert gene symbols to Entrez IDs if needed
    ids <- update_se_with_entrez_ids(data)
  } else if (ids_type == "entrez") {
    ids <- rownames(data)
  } else {
    stop("Invalid 'ids_type'. Choose either 'symbol' or 'entrez'.")
  }


  # Obtain consensus subtypes using the get.consensus.subtypes function from consensusOV
  subtypes <- get.consensus.subtypes(assay(se), ids,
                                     concordant.tumors.only = concordant.tumors.only,
                                     remove.using.cutoff = remove.using.cutoff,
                                     percentage.dataset.removed = percentage.dataset.removed)

  colData(se)$Tumor_Molecular_Subtypes <- subtypes$consensusOV.subtypes
  return(se)
}

