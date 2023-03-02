## OvRSeq: An Ovarian Cancer RNA-Seq Analysis Package

> :warning: **OvRSeq is still in development**: with great power comes great responsibility!

OvRSeq is an R package for analyzing RNA sequencing data from ovarian cancer patients. The package includes functions for quality control, normalization, differential gene expression analysis, and pathway analysis. OvaRSeq also provides options for visualization of results, such as heatmaps and volcano plots. The package is designed to be user-friendly and applicable to various types of RNA sequencing data.

### Installation

You can install the latest version of OvRSeq from GitHub with:

```
# install.packages("devtools")
devtools::install_github("icbi-lab/OvRSeq")
```

### Usage

```
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
# Quality control
qc_results <- qc(TCGA_OV)

# Normalization
norm_counts <- log2Norm(RPKMNorm(assay(TCGA_OV)))

assay(TCGA_OV) <- norm_counts


```

For more detailed information on the functions included in OvRSeq, please see the package documentation.

### License

This package is licensed under the MIT License. See the LICENSE file for details.
