## OvRSeq: An Ovarian Cancer RNA-Seq Analysis Package

OvRSeq is an R package for analyzing RNA sequencing data from ovarian cancer patients. The package includes functions for quality control, normalization, differential gene expression analysis, and pathway analysis. OvaRSeq also provides options for visualization of results, such as heatmaps and volcano plots. The package is designed to be user-friendly and applicable to various types of RNA sequencing data.

### Installation

You can install the latest version of OvRSeq from GitHub with:

```
# install.packages("devtools")
devtools::install_github("your_username/OvRSeq")
```

### Usage

```
library(OvRSeq)

# Load example data
data("TCGA_OV")

# Quality control
qc_results <- qc(TCGA_OV)

# Normalization
norm_counts <- TPMNorm(TCGA_OV, gene_lengths)

# Differential expression analysis
de_results <- diffExpr(norm_counts, col_data, "TUMOR_STAGE")

# Pathway analysis
pathway_results <- pathwayAnalysis(de_results, "hsa05200")

# Visualization
heatmap(de_results)
volcanoPlot(de_results)
```

For more detailed information on the functions included in OvRSeq, please see the package documentation.

### License

This package is licensed under the MIT License. See the LICENSE file for details.
