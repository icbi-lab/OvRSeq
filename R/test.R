#############
## Test ####
############

test <- function(){
  count_data <- read.csv("../data/TCGA-OV_RNASeq_raw.csv", row.names = 1, check.names = F)
  symbols <- rownames(count_data)
  ensembl_ids <- getEnsemblIds(symbols, org = "hsapiens_gene_ensembl")

  gene_length <- getGeneLength(ensembl_ids)

  lGene <- gene_length$external_gene_name
  count_data <- count_data[lGene,]
  countsTPM <- TPMNorm(count_data, gene_length)

  countsLog2TPM <- log2Norm(countsTPM)
}

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


icon7_dataset <- function(){
  all_data <- read.csv("../../../data/projects/2020/OvarianCancerHH/OV_R_package/Scripts/TIF_classification/Icon7_log2TPM_for_classifierTraining.csv")
  metada_sample <- all_data[c("X","TinfStatus")]
  colnames(metada_sample) <- c("Sample","TinfStatus")
  count_data <- all_data[3:27060]

  rownames(metada_sample) <- metada_sample$Sample
  rownames(count_data) <- metada_sample$Sample
  ICON7 <- SummarizedExperiment(assays=list(counts=as.matrix(t(count_data))),
                                  colData=as.data.frame(metada_sample))
  return(ICON7)

}

run_test <- function(){
  library(OvRSeq)
  se <- load_TCGA_OV()
  se <- OvRSeq(se, normalize = F)

}
