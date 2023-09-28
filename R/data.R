#' Random Forest Classifier for BRCAness
#'
#' A trained random forest classifier for predicting BRCAness status using a
#' 24-gene expression signature. This classifier was trained with multiomics
#' data from the TCGA-OV dataset and utilizes RNA sequencing (RNA-Seq) data to
#' classify patients as either having BRCAness or BRCA status.
#'
#' @format A list containing a trained random forest classifier object.
#'
#' @description
#' This data object contains a trained random forest classifier that utilizes a
#' 24-gene expression signature to classify patients into two categories: those
#' with BRCAness and those with BRCA status. The classifier was trained using
#' multiomics data, including RNA-Seq data, from the TCGA-OV dataset.
#'
#' @details
#' The random forest classifier in this data object was trained with feature
#' selection to optimize its ability to predict BRCAness status based on RNA-Seq
#' data. It is a result of a multiomics approach aimed at accurately classifying
#' patients.
#'
#' @usage data(brcaness_classifier)
#'
#' @examples
#' # Load the trained BRCAness classifier
#' data(brcaness_classifier)
#'
#' # Predict BRCAness status for a new patient using the classifier
#' new_patient_data <- ...  # Replace with new patient's RNA-Seq data
#' prediction <- predict(brcaness_classifier, new_patient_data)
#'
#' @docType data
#' @name brcaness_classifier
NULL

#' TCGA-OV RNA-Seq Dataset
#'
#' A SummarizedExperiment object containing raw RNA-Seq counts from the TCGA-OV dataset,
#' along with metadata columns including AGE, TUMOR_GRADE, CLINICAL_STAGE,
#' TUMOR_RESIDUAL_DISEASE, and BRCAness.
#'
#' @format
#' A SummarizedExperiment object.
#'
#' @description
#' This dataset contains raw RNA-Seq counts from the TCGA-OV dataset. It includes
#'  metadata columns that provide additional information about the samples,
#'  including patient age, tumor grade, clinical stage, tumor residual disease
#'  status, and BRCAness classification.
#'
#' @details
#' The TCGA-OV dataset is stored as a SummarizedExperiment object, which provides
#' an organized and efficient structure for working with high-throughput genomics data.
#' The metadata columns allow for the exploration of clinical and molecular attributes
#' associated with the samples.
#'
#' @seealso
#' \code{SummarizedExperiment} for more information on SummarizedExperiment objects.
#'
#' @usage data(TCGA_OV)
#'
#' @examples
#' # Load the TCGA-OV dataset
#' data(TCGA_OV)
#'
#' # Access metadata columns
#' metadata(TCGA_OV)$AGE
#' metadata(TCGA_OV)$TUMOR_GRADE
#'
#' @docType data
#' @name TCGA_OV
NULL

#' BRCAness Gene Signature
#'
#' A gene signature developed through feature selection to optimally classify
#' TCGA-OV samples into BRCAness and non-BRCAness categories. The signature
#' consists of gene expression values (2TPM+1 normalized) for 24 genes. These
#' genes were selected based on their importance in three different machine
#' learning models (Random Forest, Ada Boost, and Gradient Boosting) after
#' recursive feature elimination. Only genes that were among the top 50 most
#' important features in at least two of the three models were included in
#' this signature.
#'
#' @format
#' A vector containing the gene symbol for 24 genes.
#'
#' @description
#' The BRCAness gene signature is a result of an extensive feature selection
#' process using machine learning models. It represents a subset of genes that
#' are critical for discriminating between BRCAness and non-BRCAness samples
#' based on gene expression data from the TCGA-OV dataset.
#'
#' @seealso
#' \code{TCGA_OV} for the dataset used for feature selection and training.
#'
#' @usage data(brcaness_signature)
#'
#' @examples
#' # Load the BRCAness gene signature
#' data(brcaness_signature)
#'
#' # Access gene expression values for the first gene
#' brcaness_signature
#'
#' @docType data
#' @name brcaness_signature
NULL

#' Immune Signatures
#'
#' A GMT (Gene Matrix Transposed) file containing well-defined immune-related
#' gene signatures. Each signature consists of a list with the name of the
#' signature and a vector of gene symbols representing the genes that constitute
#' the signature. These signatures are used to characterize immune-related
#' processes, including T cell inflammation, IFN gamma signature, cytolytic
#' activity, and cytotoxic T lymphocyte function.
#'
#' @format
#' A list of lists, where each sublist contains:
#'   - `Name`: The name of the immune signature.
#'   - `Genes`: A character vector of gene symbols representing the signature genes.
#'
#' @description
#' The `immune_signatures` object contains well-defined immune-related gene
#' signatures in GMT format. These signatures are valuable for characterizing
#' immune processes and activities in biological data, including gene expression
#' profiles.
#'
#' @seealso
#' For more information on GMT files and gene set enrichment analysis (GSEA),
#' refer to the relevant literature and software documentation.
#'
#' @usage data(immune_signatures)
#'
#' @examples
#' # Load the immune signatures
#' data(immune_signatures)
#'
#' # Access the genes in the "T cell inflammation" signature
#' t_cell_inflammation_genes <- immune_signatures$T_cell_inflammation$Genes
#'
#' @docType data
#' @name immune_signatures
NULL

#' Tumor Immune Phenotype Signature
#'
#' A gene signature developed based on digital pathology and transcriptome analysis
#' to classify tumor immune phenotypes (infiltrate, excluded, desert) in ovarian cancer.
#' This signature was derived from a classification model using 157 genes that describe
#' the presence and position of CD8+ T cells relative to the tumor center or margin
#' in the ICON7 cohort.
#'
#' @format
#' A list containing gene symbols representing the 148 genes in the tumor immune phenotype signature.
#'
#' @description
#' The `tumor_immune_phenotype_signature` object contains a gene signature that classifies
#' ovarian cancer tumor immune phenotypes. It was developed using integrated digital pathology
#' and transcriptome analysis, with a focus on CD8+ T cell presence and position within the tumor.
#'
#' @references
#' Desbois M, Udyavar AR, Ryner L, Kozlowski C, Guan Y, Dürrbaum M, et al. Integrated digital
#' pathology and transcriptome analysis identifies molecular mediators of T-cell exclusion in
#' ovarian cancer. Nat Commun. 2020;11:5583.
#'
#' @usage data(tumor_immune_phenotype_signature)
#'
#' @examples
#' # Load the tumor immune phenotype signature
#' data(tumor_immune_phenotype_signature)
#'
#' # Access the gene symbols in the signature
#' signature_genes <- tumor_immune_phenotype_signature
#'
#' @docType data
#' @name tumor_immune_phenotype_signature
NULL

