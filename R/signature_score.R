

#' Calculate enrichment scores for immune signatures
#'
#' This function calculates enrichment scores for a list of immune signatures
#' using the GSVA method. The user can provide their own gene list or use one
#' of the pre-defined ones. The function returns an updated SummarizedExperiment
#' object with the enrichment scores added to the colData.
#'
#' @param se A SummarizedExperiment object containing RNA sequencing data
#' @param method A character string indicating the method to use for calculating enrichment scores. Default is "ssgsea".
#' @param genelist A character vector with the gene list for the immune signature.
#'
#' @return A SummarizedExperiment object with the enrichment scores added to the colData.
#' @importFrom GSVA gsva
#' @export
immune_signature_score <- function(se, method = "ssgsea", genelist ){
  count_data <- assay(se)
  results <- GSVA::gsva(count_data,
                 genelist,
                 method = 'ssgsea',
                 ssgsea.norm=FALSE)
  results <- t(results)
  colData(se) <- cbind(colData(se), results)
  return(se)
}

