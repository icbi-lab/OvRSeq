#' Perform Single-Sample Gene Set Enrichment Analysis (ssGSEA) on a SummarizedExperiment
#'
#' This function performs ssGSEA on a given `SummarizedExperiment` object using custom gene sets provided in GMT format.
#' It leverages the `GSVA` package to compute enrichment scores for each sample in the experiment.
#' GSVA is a non-parametric, unsupervised method for estimating variation of gene set enrichment
#' through the samples of an expression dataset.
#'
#' @param se A `SummarizedExperiment` object containing the expression data.
#' @param gmt_list A list where each element is a character vector of gene symbols representing a gene set.
#'        The names of the list elements are used as gene set names.
#' @return A matrix of enrichment scores with gene sets as rows and samples as columns.
#'
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object and 'gmt_list' is your list of gene sets
#' enrichment_scores <- ssGSEA_OV_custom(se, gmt_list)
#'
#' @importFrom GSVA gsva
#' @export
#' @references
#' Hänzelmann S, Castelo R, Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.”
#' BMC Bioinformatics, 14, 7. <doi:10.1186/1471-2105-14-7>
#' For a full list of citations, use 'citation("GSVA")' in R.
ssGSEA_OV_custom <- function(se, gmt_list) {
  # Perform SS-GSEA using the GSVA package
  # Notice: the user should ensure that the gene identifiers in the expression data and gene sets match
  results <- gsva(se, gmt_list, method = 'ssgsea', ssgsea.norm = FALSE)
  colData(se) <- cbind(colData(se), t(assay(results)))
  # Return the matrix of enrichment scores
  return(se)
}
