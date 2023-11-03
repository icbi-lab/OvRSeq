

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

#' Compute average expression values for gene signatures and enrich colData of a SummarizedExperiment object
#'
#' Given a `SummarizedExperiment` object and a matrix or data frame of gene signatures, 
#' this function computes the average expression values for each gene signature and enriches 
#' the `colData` of the `SummarizedExperiment` object with the computed values.
#'
#' @param se A `SummarizedExperiment` object with gene expression data.
#' @param gene_sig A matrix or data frame where each column represents a gene signature and each cell contains a gene symbol.
#'        Empty cells should be represented as empty strings. 
#' @param genesetname A character string that serves as a prefix for the new columns added to the `colData` of the `SummarizedExperiment` object.
#' 
#' @return A `SummarizedExperiment` object enriched with new columns in `colData`, each representing 
#'         the average expression value of a given gene signature for each sample.
#' 
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object and 'gene_sig' is your gene signatures matrix
#' enriched_se <- avg_expression_for_signature_se(se, gene_sig, "MyGeneSet")
#'
#' @export
avg_expression_for_signature_se <- function(se, gene_sig, genesetname) {
  
  # Check if the genesetname is already present in colData, if so provide a warning
  if (genesetname %in% colnames(colData(se))) {
    warning("The genesetname is already present in colData. Results will overwrite existing values.")
  }
  
  # Extract signatures from gene_sig
  sigs <- colnames(gene_sig)
  
  # Extract expression data from the SummarizedExperiment object
  expression <- assay(se)
  
  # For each signature, compute average expression
  for (i in seq_along(sigs)) {
    genlist <- gene_sig[, i]
    
    # Filter out empty entries
    genlist <- genlist[genlist != ""]
    
    # Match rownames from se with genlist to handle any missing genes
    matched_genes <- rownames(se)[rownames(se) %in% genlist]
    
    # Get expression values for genes in the signature
    expression_genes <- expression[matched_genes, , drop = FALSE]
    
    # Compute average expression across genes for each sample
    av_expression <- colMeans(expression_genes, na.rm = TRUE)
    
    # Add this to colData of the SummarizedExperiment
    colData(se)[[paste0(genesetname, "_", sigs[i])]] <- av_expression
  }
  
  return(se)
}

