#' Deconvolute Immune Cell Fractions from Gene Expression Data
#'
#' This function deconvolutes immune cell fractions from gene expression data using methods from the `immunedeconv` package.
#'
#' @param se A Summarized Experiment of gene expression values (logTPM+1), where rows are genes and columns are samples.
#' @param method A character string indicating the deconvolution method to be used. See details for available methods.
#' @return A list with deconvolution results.
#' @importFrom immunedeconv deconvolute
#' @examples
#' # data <- ... # Your gene expression data (logTPM+1)
#' # res <- deconvolute_immune(data, method="timer")
#' @export
#' @details
#' Available methods and their respective licensing and citations are as follows:
#' 
#' \tabular{lll}{
#' method \tab license \tab citation \cr
#' quanTIseq \tab free (BSD) \tab Finotello et al. (2019) Genome Medicine \cr
#' TIMER \tab free (GPL 2.0) \tab Li et al. (2016) Genome Biology \cr
#' CIBERSORT \tab free for non-commerical use only \tab Newman et al. (2015) Nature Methods \cr
#' MCPCounter \tab free (GPL 3.0) \tab Becht et al. (2016) Genome Biology \cr
#' xCell \tab free (GPL 3.0) \tab Aran et al. (2017) Genome Biology \cr
#' EPIC \tab free for non-commercial use only \tab Racle et al. (2017) ELife \cr
#' }
#' 
#' Note: While `immunedeconv` itself is free (BSD), you may need to obtain a license to use individual methods.
#' Please ensure you cite both this package and the method(s) you use in your work.
#'
#' @references
#' Sturm, G., et al. (2019) Bioinformatics. \url{https://doi.org/10.1093/bioinformatics/btz363}
#' 
deconvolute_immune <- function(se, method) {
  indications <- rep('OV', times=ncol(data))
  data <- assay(se)
  if (method == "timer") {
    res <- deconvolute(data, method=method, indications=indications)
  } else {
    res <- deconvolute(data, method=method)
  }
   lCell <- unlist(lapply(res$cell_type,function(x){paste0(x,"|", method)}))
   res$cell_type <- NULL
   res <- t(res)
   colnames(res) <- lCell
   colData(se) <- cbind(colData(se), res)
   return(se)
}
