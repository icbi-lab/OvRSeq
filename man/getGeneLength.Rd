% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normalization.R
\name{getGeneLength}
\alias{getGeneLength}
\title{Retrieve gene length from Ensembl IDs}
\usage{
getGeneLength(ensembl_ids, org = "hsa")
}
\arguments{
\item{ensembl_ids}{A data frame containing two columns: \code{ensembl_gene_id}
and \code{external_gene_name}. The former contains the Ensembl gene IDs of the genes
for which the gene length should be retrieved, while the latter contains their
corresponding external gene symbols.}

\item{org}{A character string indicating the organism for which gene length should be
retrieved. The default value is "hsa", corresponding to the human genome.}
}
\value{
A data frame containing three columns: \code{ensembl_gene_id},
\code{external_gene_name}, and \code{length}. The length column contains the length of
each gene in base pairs (bp).
}
\description{
This function retrieves the length of genes based on their Ensembl IDs. The gene length
is defined as the sum of the lengths of all exons of the gene. The function uses the
\code{exonsBy} function from the \pkg{ensembldb} package to obtain exon information,
and the \pkg{EnsDb.Hsapiens.v86} package for the human reference database.
}
