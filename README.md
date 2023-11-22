# OvRSeq: An Ovarian Cancer RNA-Seq Analysis Package

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/doseRider)](https://cran.r-project.org/package=ovrseq)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
> :warning: **OvRSeq is still in development**: with great power comes great responsibility!

Cancer immunotherapy has shown limited efficacy in high-grade serous ovarian cancer (HGSOC), prompting clinical trials combining immunotherapy with poly-ADP-ribose polymerase (PARP) inhibitors. Identifying predictive biomarkers is crucial. We introduce a novel diagnostic algorithm and an R package, OvRSeq, for comprehensive RNA sequencing-based characterization of HGSOC samples. Leveraging a 24-gene expression signature, OvRSeq classifies samples for BRCAness, a homologous recombination repair deficiency (HRD) phenotype. It thoroughly profiles the immune environment, including molecular subtype, tumor-immune phenotype, immunophenoscore (IPS), established and novel immune signatures (e.g., IFNg gene set, T cell inflammation, exhaustion, cytolytic activity, immune checkpoints), and estimates infiltrated immune cells via deconvolution methods. OvRSeq also provides an indicator for combination therapy response by assessing BRCAness within an immune environment balanced between activation and suppression.

## Dependencies

OvRseq has dependencies on packages outside CRAN. Please make sure you install them before installing it.


```
remotes::install_github("GfellerLab/EPIC")
BiocManager::install("consensusOV")
```

## Installation

You can install the latest version of OvRSeq from GitHub with:

```
# install.packages("devtools")
devtools::install_github("icbi-lab/OvRSeq")
```
Certainly! Below is an example of how you might rewrite the README.md for the OvRSeq package to include descriptions of the implemented functions:


## Features and Functions

Key features provided by OvRSeq include:

- Evaluation of BRCAness status by leveraging a predictive model developed from the TCGA-OV dataset, which helps in identifying the likelihood of BRCAness in ovarian cancer samples.

- Determination of tumor immune phenotypes by assessing the spatial distribution patterns of CD8+ T cells within the tumor microenvironment, aiding in the understanding of immune evasion mechanisms in cancer.

- Molecular subtype classification through a consensus approach that categorizes ovarian cancer into distinct subtypes, enabling tailored therapeutic strategies.

- Immune profile characterization with the calculation of the Immunophenoscore (IPS), which includes assessments of various immune cell types and their activation states.

- Immune deconvolution from gene expression data, which parses out the proportions of different immune cell types, providing insights into the immune landscape of ovarian cancer samples.

- Calculation of enrichment scores for multiple immune-related gene signatures, employing the GSVA method to determine the relative activation of different immune pathways.

- Averaging of expression levels across smaller, defined immune-related gene sets, offering a precise quantification of signature expressions in a given dataset.

OvRSeq effectively serves as a valuable resource for researchers focusing on the molecular and immunological characterization of ovarian cancer, supporting the drive towards precision medicine by furnishing a multi-dimensional analysis of gene expression data.


## Example Usage

Here's how to calculate the average expression for a signature and enrich `colData`:

```r
library(OvRSeq)
#Load TCGA-OV dataset 
load_TCGA_OV()
```

To perform the main analysis :

```r
se <- OvRSeq(TCGA_OV, normalize = F)
```
Plot results
```{r}
plot_ggmarginal(se, x_var = "CD8 Jiang" , y_var =  "Monocyte|quantiseq" , color_var = "BRCAness")
```

For a full list of functions and their descriptions, consult the package documentation.

## Development and Contributions

OvRSeq is under active development. We welcome contributions and suggestions from the community. Please feel free to submit issues and pull requests through our [GitHub repository](https://github.com/icbi-lab/OvRSeq).

## License

OvRSeq is released under the [MIT License](https://opensource.org/licenses/MIT).


