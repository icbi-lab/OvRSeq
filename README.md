## OvRSeq: An Ovarian Cancer RNA-Seq Analysis Package

> :warning: **OvRSeq is still in development**: with great power comes great responsibility!

Cancer immunotherapy has shown limited therapeutic efficacy in high-grade serous ovarian cancer (HGSOC). To overcome these limitations, a number of clinical trials combining immunotherapies with poly-ADP-ribose polymerase (PARP) inhibitors are currently underway. Hence there is an urgent need to identify predictive biomarkers. Here we propose a novel diagnostic algorithm and provide a companion R package (**OvRSeq**) for comprehensive RNA sequencing based characterization of HGSOC samples. Based on a 24 gene expression signature samples will be classified according to BRCAness, which is an homologous recombination repair deficiency (HRD) phenotype beyond BRCA1/2 mutations. The immune environment will be comprehensively characterized by determination of the molecular subtype (IMR, DIF, PRO, MES), the tumor-immune phenotype (Infilterated, Excluded, Desert), the immunophenoscore (IPS), enrichment of a number of well established and novel immune signatures, such as IFNg gene set, T cell inflammation, T cell exhaustion, cytolytic activity, immune checkpoints as well as estimation of infiltrated immune cells using various deconvolution methods. Further an indicator for response to combination therapy that indicates BRCAness with a favorable balance between activated and suppressive immune environment is provided.

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
