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
#' data(BRCAness_signature)
#'
#' # Train a random forest model using the BRCAness label
#' rf_model <- train_rf(TCGA_OV, "BRCAness", rownames(BRCAness_signature))
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


#' Load pre-trained BRCAness random forest model
#'
#' This function loads a pre-trained random forest model for BRCAness classification from a saved .rda file.
#'
#' @return A trained random forest model object
#' @examples
#' brcaness_model <- load_brcaness_classifier()
load_brcaness_classifier <- function(){
  load("../data/rf_model.rda")
  return(rf_model)
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
classify_brcaness <- function(se, brcaness_classifier=NA, brcaness_signature=NA) {

  if (is.na(brcaness_classifier)){
    print("No BRCAness classifier has been provieded. Using pre-trained RF model")
    brcaness_classifier <- load_brcaness_classifier()
  }

  if (is.na(brcaness_signature)){
    print("No BRCAness classifier has been provieded. Using pre-selected signature")
    brcaness_signature <- scan("../data/BRCAnessSignature.csv", what = character())
    }

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

