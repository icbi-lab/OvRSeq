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

## Usage

```{r}
#  Load the OvRSeq package
library(OvRSeq)

# Load example data
data("TCGA_OV")
data("brcaness_signature")

# Classify BRCAness 
TCGA_OV <- OvRSeq::classify_brcaness(TCGA_OV, rf_model, brcaness_signature)

# Check results
table(colData(TCGA_OV)$BRCAness)

# BRCAness noBRCAness 
#      107        119 
=======

# Obtain immune signature score
## Load the geneset
data("immune_signatures")
TCGA_OV <- immune_signature_score(TCGA_OV, immune_signatures)
```

For more detailed information on the functions included in OvRSeq, please see the package documentation.

## Contributing

Contributions to `OvRSeq` are welcome. If you find any bugs, have feature requests, or want to contribute improvements or new features, please open an issue or submit a pull request on GitHub.

## License

This package is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
