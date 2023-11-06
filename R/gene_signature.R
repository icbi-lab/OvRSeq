#' Load BRCAness Gene Signature
#'
#' Loads the BRCAness gene signature from the "brcaness_signature.rda" data file
#' and returns it as a vector with gene symbol for the for 24 genes.
#'
#' @return
#' A vector with gene symbol for the BRCAness gene signature.
#'
#' @usage
#' load_brcaness_signature()
#'
#' @examples
#' # Load the BRCAness gene signature
#' brcaness_signature <- load_brcaness_signature()
#'
#' # Access gene expression values for the first gene
#' brcaness_signature$Gene1
#'
#' @export
load_brcaness_signature <- function() {
  load("data/brcaness_signature.rda")
  if (exists("brcaness_signature", inherits = FALSE)) {
    return(brcaness_signature)
  } else {
    stop("brcaness_signature not found in the loaded data.")
  }
}

#' Load Immune Signatures
#'
#' Loads the immune signatures from the "immune_signatures.rda" data file and
#' returns them as a list of lists, where each sublist contains the name of
#' an immune signature and a vector of gene symbols representing the genes that
#' constitute the signature.
#'
#' @return
#' A list of lists, where each sublist contains:
#'   - `Name`: The name of the immune signature.
#'   - `Genes`: A character vector of gene symbols representing the signature genes.
#'
#' @usage
#' load_immune_signatures()
#'
#' @examples
#' # Load the immune signatures
#' immune_signatures <- load_immune_signatures()
#'
#' # Access the genes in the "T cell inflammation" signature
#' t_cell_inflammation_genes <- immune_signatures$T_cell_inflammation$Genes
#'
#' @export
load_immune_signatures <- function() {
  load("data/immune_signatures.rda")
  if (exists("immune_signatures", inherits = FALSE)) {
    return(immune_signatures)
  } else {
    stop("immune_signatures not found in the loaded data.")
  }
}

#' Load Small Immune Signatures
#'
#' Loads the small immune signatures from the "small_immune_signatures.rda" data file and
#' returns them as a list of lists. Each sublist corresponds to an immune signature that
#' contains fewer than 10 genes, making it suitable for analysis methods like z-score computations
#' that may be more appropriate for smaller gene sets than single-sample gene set variance analysis.
#'
#' @return
#' A list of lists, where each sublist contains:
#'   - `Name`: The name of the small immune signature.
#'   - `Genes`: A character vector of gene symbols representing the genes in the signature.
#'
#' @usage
#' small_immune_signatures <- load_small_immune_signatures()
#'
#' @examples
#' # Load the small immune signatures
#' small_immune_signatures <- load_small_immune_signatures()
#'
#' # Access the genes in the "T cell inflammation" signature
#' ifng_ayers <- small_immune_signatures$'IFNG Ayers'
#'
#' @export
load_small_immune_signatures <- function() {
  load("data/small_immune_signatures.rda")
  if (exists("small_immune_signatures", inherits = FALSE)) {
    return(small_immune_signatures)
  } else {
    stop("small_immune_signatures not found in the loaded data.")
  }
}
