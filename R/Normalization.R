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


#' Get Ensembl IDs from Gene Symbols
#'
#' Takes a vector of gene symbols and returns a data frame of Ensembl gene IDs
#' and gene symbols using the biomaRt package.
#'
#' @param gene_symbols A character vector of gene symbols to look up
#' @param org A string specifying the organism (default: "hsapiens_gene_ensembl")
#' @return A data frame with two columns: "ensembl_gene_id" and "external_gene_name"
#' @importFrom biomaRt useMart getBM
#' @examples
#' getEnsemblIds(c("BRCA1", "TP53"))
#' @export
getEnsemblIds <- function(gene_symbols, org = "hsapiens_gene_ensembl") {
  mart <- useMart("ensembl", dataset = org)
  attributes <- c("ensembl_gene_id", "external_gene_name")
  filters <- "external_gene_name"
  ensembl_ids <- getBM(attributes = attributes, filters = filters,
                       values = gene_symbols, mart = mart)
  return(ensembl_ids)
}



#' Retrieve gene length from Ensembl IDs
#'
#' This function retrieves the length of genes based on their Ensembl IDs. The gene length
#' is defined as the sum of the lengths of all exons of the gene. The function uses the
#' \code{exonsBy} function from the \pkg{ensembldb} package to obtain exon information,
#' and the \pkg{EnsDb.Hsapiens.v86} package for the human reference database.
#'
#' @param ensembl_ids A data frame containing two columns: \code{ensembl_gene_id}
#' and \code{external_gene_name}. The former contains the Ensembl gene IDs of the genes
#' for which the gene length should be retrieved, while the latter contains their
#' corresponding external gene symbols.
#' @param org A character string indicating the organism for which gene length should be
#' retrieved. The default value is "hsa", corresponding to the human genome.
#'
#' @return A data frame containing three columns: \code{ensembl_gene_id},
#' \code{external_gene_name}, and \code{length}. The length column contains the length of
#' each gene in base pairs (bp).
#'
#' @importFrom ensembldb exonsBy
#' @importFrom IRanges reduce
#' @import EnsDb.Hsapiens.v86
#' @export
getGeneLength <- function(ensembl_ids, org = "hsa"){

  exons = exonsBy(EnsDb.Hsapiens.v86, by="gene")
  exons = reduce(exons)
  len = sum(width(exons))

  insect = intersect(ensembl_ids$ensembl_gene_id,names(len))
  geneLengths = len[insect]
  rownames(ensembl_ids) = ensembl_ids$ensembl_gene_id
  ensembl_ids$length = NA
  ensembl_ids[insect,]$length = as.vector(unlist(geneLengths))
  ensembl_ids$KB <- ensembl_ids$length/1000
  ensembl_ids = ensembl_ids[!is.na(ensembl_ids$length),]
  return(ensembl_ids)
}

#' TPM normalization of count data
#'
#' Takes a count data matrix and a gene length data frame, and applies TPM normalization.
#'
#' @param count_data A matrix of count data with genes in rows and samples in columns.
#' @param gene_length A vector with gene IDs as names and transcript lengths in kilobases (KB).
#'
#' @return A TPM-normalized matrix with the same dimensions as \code{count_data}.
#'
#' @examples
#' # Load example data from SummarizedExperiment package
#' data("example_se")
#' # Extract count data and gene length data
#' counts <- assay(example_se)
#' gene_length <- rowData(example_se)$gene_length
#' # Apply TPM normalization
#' normalized_counts <- TPMNorm(counts, gene_length)
#'
#' @export
TPMNorm <- function(count_data,gene_length){
  # Create new hnsc object, don't overwrite the previous
  countsTPM <- count_data

  # Divide each gene by transcript length
  countsTPM <- apply( countsTPM, 2, function(x){ x / gene_length } )

  # Divide by the transcript length
  countsTPM <- apply( countsTPM, 2, function(x) { x / sum(x) * 1E6} )

  rownames(countsTPM) <- rownames(count_data)
  return(countsTPM)
}


#' RPKM normalization of count data
#'
#' Takes a count data matrix, and applies RPKM normalization.
#'
#' @param count_data A matrix of count data with genes in rows and samples in columns.
#'
#' @return A RPKM-normalized matrix with the same dimensions as \code{count_data}.
#'
#' @examples
#' # Load example data from SummarizedExperiment package
#' data("example_se")
#' # Extract count data and gene length data
#' counts <- assay(example_se)
#' # Apply RPKM normalization
#' normalized_counts <- RPKMNorm(counts)
#'
#' @export
RPKMNorm <- function(count_data){
  countPKM <- apply( count_data, 2, function(x) { x / sum(x) * 1E6} )
  return(countPKM)
}

#' Log2 normalization of count data
#'
#' This function performs log2 normalization of count data. It adds 1 to each count value and takes the logarithm base 2 of the resulting value.
#'
#' @param count_data A matrix of count data with genes in rows and samples in columns..
#'
#' @return A matrix of count data with genes in rows and samples in columns, containing the log2 normalized count data.
#'
#' @examples
#' data("example_counts")
#' data("example_gene_length")
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = example_counts))
#' log2_norm_data <- log2Norm(se)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
log2Norm <- function(count_data){
  # TPM
  count_data_log2 <- count_data
  count_data_log2 <- log(count_data_log2 + 1, 2)
  return(count_data_log2)
}


