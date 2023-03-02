
test_data <- function(){
  all_data <- read.csv("../../data/Matrix_log2tpm.csv", check.names = F, row.names = 1)
  lCol_data <- c("AGE","TUMOR_GRADE","CLINICAL_STAGE","TUMOR_RESIDUAL_DISEASE","BRCAness")

  OV <- read.csv("../data/TCGA-OV_RNASeq_raw.csv", check.names = F, row.names = 1)
  lGene <- intersect(rownames(OV), colnames(all_data))
  count_data_log_tpm <- all_data[lGene]
  col_data <- all_data[lCol_data]
  count_data_log_tpm = t(count_data_log_tpm)
  colnames(count_data_log_tpm) = rownames(all_data)

  library(SummarizedExperiment)
  #Create SummarizedExperiment object

  tcga_ov <- SummarizedExperiment(assays=list(counts=as.matrix(count_data_log_tpm)),
                             colData=as.data.frame(col_data))


  save(se, file = "../data/TCGA_OV.rda")
}
