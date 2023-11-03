# OvRSeq: An Ovarian Cancer RNA-Seq Analysis Package

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/doseRider)](https://cran.r-project.org/package=ovrseq)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
> :warning: **OvRSeq is still in development**: with great power comes great responsibility!

Cancer immunotherapy has shown limited efficacy in high-grade serous ovarian cancer (HGSOC), prompting clinical trials combining immunotherapy with poly-ADP-ribose polymerase (PARP) inhibitors. Identifying predictive biomarkers is crucial. We introduce a novel diagnostic algorithm and an R package, OvRSeq, for comprehensive RNA sequencing-based characterization of HGSOC samples. Leveraging a 24-gene expression signature, OvRSeq classifies samples for BRCAness, a homologous recombination repair deficiency (HRD) phenotype. It thoroughly profiles the immune environment, including molecular subtype, tumor-immune phenotype, immunophenoscore (IPS), established and novel immune signatures (e.g., IFNg gene set, T cell inflammation, exhaustion, cytolytic activity, immune checkpoints), and estimates infiltrated immune cells via deconvolution methods. OvRSeq also provides an indicator for combination therapy response by assessing BRCAness within an immune environment balanced between activation and suppression.

## Installation

You can install the latest version of OvRSeq from GitHub with:

```
# install.packages("devtools")
devtools::install_github("icbi-lab/OvRSeq")
```
Certainly! Below is an example of how you might rewrite the README.md for the OvRSeq package to include descriptions of the implemented functions:


## Features and Functions

OvRSeq includes the following key functions:

- `avg_expression_for_signature_se`: Computes average gene expression values and enriches `colData` in `SummarizedExperiment` objects.

- `brcaness_classifier`: Contains a trained Random Forest model for predicting BRCAness based on a 24-gene expression signature derived from the TCGA-OV dataset.

- `calculateIPS`: Calculates the Immunophenoscore (IPS) and its component scores from a `SummarizedExperiment` object.

- `classify_brcaness`: Uses the BRCAness classifier to categorize samples based on their likelihood of BRCAness.

- `deconvolute_immune`: Deconvolutes immune cell fractions from gene expression data using the `immunedeconv` package.

- `get_consensus_ov_subtypes`: Implements a consensus classifier to determine molecular subtypes of high-grade serous ovarian cancer.

- `immune_signature_score`: Calculates enrichment scores for various immune signatures using the gene set variation analysis (GSVA) method.

- `load_TCGA_OV`: Provides an interface to load the TCGA-OV dataset as a `SummarizedExperiment` object.

- `ssGSEA_OV_custom`: Performs single-sample Gene Set Enrichment Analysis (ssGSEA) using custom gene sets and the `GSVA` package.

- `tumor_immune_phenotype_signature`: A gene signature for classifying tumor immune phenotypes in ovarian cancer based on the spatial distribution of CD8+ T cells.

## Example Usage

Here's how to calculate the average expression for a signature and enrich `colData`:

```r
library(OvRSeq)
# Assuming 'se' is a SummarizedExperiment object
se <- avg_expression_for_signature_se(se)
```

To classify samples for BRCAness:

```r
brcaness_status <- classify_brcaness(se)
```

For a full list of functions and their descriptions, consult the package documentation.

## Development and Contributions

OvRSeq is under active development. We welcome contributions and suggestions from the community. Please feel free to submit issues and pull requests through our [GitHub repository](https://github.com/icbi-lab/OvRSeq).

## License

OvRSeq is released under the [MIT License](https://opensource.org/licenses/MIT).


