#' Update Gene Symbols to Entrez IDs in a SummarizedExperiment
#'
#' This function takes a SummarizedExperiment object with gene symbols and updates the row names with corresponding Entrez IDs.
#'
#' @param se A SummarizedExperiment object.
#' @return A SummarizedExperiment object with Entrez IDs as row names.
#' @import org.Hs.eg.db
#' @import SummarizedExperiment
#' @importFrom AnnotationDbi mapIds
#' @examples
#' # Assuming you have a SummarizedExperiment object named 'se'
#' se_updated <- update_se_with_entrez_ids(se)
#' rownames(se_updated)
update_se_with_entrez_ids <- function(se) {
  # Extract gene symbols from the SummarizedExperiment
  gene_symbols <- rownames(se)

  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = gene_symbols,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

  # Update row names of the SummarizedExperiment with Entrez IDs

  return(entrez_ids)
}

