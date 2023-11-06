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

#' Small Immune Signatures
#'
#' A GMT (Gene Matrix Transposed) file containing concise immune-related
#' gene signatures. Each signature is limited to fewer than 10 genes, which
#' is why they are termed "small." These limited gene sets can provide more
#' precise z-score calculations for immune-related processes as opposed to
#' single-sample gene set enrichment variances.
#'
#' Each signature within the file facilitates the characterization of
#' specific immune-related activities, such as T cell inflammation, interferon-gamma
#' response, cytolytic activity, cytotoxic T lymphocyte response, and various
#' components of the immune system. Splitting larger signatures into smaller,
#' more manageable sets can yield a clearer understanding of the biological
#' data by focusing on the most relevant and influential genes.
#'
#' @format
#' A list of lists, where each sublist contains:
#'   - `Name`: The name of the immune signature.
#'   - `Genes`: A character vector of gene symbols, each representing fewer than
#'     10 signature genes for focused analysis.
#'
#' @description
#' The `small_immune_signatures` object contains well-curated, small immune-related
#' gene signatures in GMT format, suitable for refined gene expression analysis.
#' These small signatures are derived from larger gene sets, trimmed for optimal
#' focus and precision in immune process characterization.
#'
#' @seealso
#' For more detailed information on the use of small gene sets for z-score calculations
#' in immune signature analysis, as well as general GMT file structure and gene set
#' enrichment analysis (GSEA), refer to the relevant literature and software documentation.
#'
#' @usage data(small_immune_signatures)
#'
#' @examples
#' # Load the small immune signatures
#' data(small_immune_signatures)
#'
#' # Access the genes in the "T cell inflammation" signature
#' t_cell_inflammation_genes <- small_immune_signatures$T_cell_inflammation$Genes
#'
#' @docType data
#' @name small_immune_signatures
#' @keywords datasets
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

#' Tumor Infiltration Status Classifier
#'
#' A Random Forest classifier trained on gene expression data from the ICON7 trial.
#' This classifier is designed to classify the tumor immune infiltration status
#' into categories such as 'infiltrated', 'excluded', or 'desert' based on a
#' gene signature. The model is trained using `randomForest` package and can be
#' used to predict the infiltration status of ovarian cancer samples.
#'
#' @format An object of class `randomForest` (inherits from `list`), representing
#'   a fitted Random Forest model.
#'
#' @details
#' The classifier was trained using a Random Forest algorithm with 300 trees
#' on the ICON7 dataset, which contains gene expression data for ovarian cancer
#' samples. The gene signature used for the training consists of genes identified
#' as common between the SummarizedExperiment datasets and the tumor immune phenotype
#' signature.
#'
#' The `classifier_infiltration_status` is saved as an R object of class `randomForest`,
#' which contains the entire fitted model object. This model can be used to predict
#' infiltration status in other datasets by supplying the appropriate gene expression
#' data for the genes included in the signature.
#'
#' The classifier should be applied only to data that has been preprocessed in the same
#' manner as the ICON7 trial data to ensure the validity of the predictions.
#'
#' @usage
#' data(classifier_infiltration_status)
#'
#' @references
#' Desbois M, Udyavar AR, Ryner L, Kozlowski C, Guan Y, Dürrbaum M, et al. Integrated digital
#' pathology and transcriptome analysis identifies molecular mediators of T-cell exclusion in
#' ovarian cancer. Nat Commun. 2020;11:5583.
#'
#' @examples
#' data(classifier_infiltration_status)
#' # Suppose `new_data` is a matrix of gene expression values with rows as genes and columns as samples:
#' predicted_status <- predict(classifier_infiltration_status, new_data)
#' @name classifier_infiltration_status
NULL

#' Immunophenoscore (IPS) Genes and Weights
#'
#' A dataset containing genes and corresponding weights used to calculate the
#' Immunophenoscore (IPS), which is a measure of the immune landscape of a tumor.
#'
#' @format
#' A data frame with the following columns:
#' \describe{
#'   \item{GENE}{(character) Official gene symbol.}
#'   \item{WEIGHT}{(numeric) Weight of the gene in the IPS calculation.}
#'   \item{CATEGORY}{(character) The immunological category to which the gene belongs
#'   (e.g., MHC, CP, EC, SC, MDSC, TREG, AZ).}
#' }
#'
#' @source
#' Immunophenogram project: https://github.com/icbi-lab/Immunophenogram
#'
#' @references
#' Charoentong P, Finotello F, Angelova M, Mayer C, Efremova M, Rieder D, Hackl H,
#' Trajanoski Z. Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype
#' Relationships and Predictors of Response to Checkpoint Blockade. Cell Rep.
#' 2017 Jan 3;18(1):248-262. doi: 10.1016/j.celrep.2016.12.019. PMID: 28052254.
#'
#' @usage
#' data(IPS_genes)
#'
#' @examples
#' # Load the IPS genes and weights
#' data(IPS_genes)
#'
#' # Access the gene symbols and weights
#' ips_genes <- IPS_genes$GENE
#' ips_weights <- IPS_genes$WEIGHT
#'
#' @docType data
#' @name IPS_genes
NULL
