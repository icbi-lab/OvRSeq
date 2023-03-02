

# run single sample gensetenrichment analysis on built in genesets
ssGSEA_OV <- function(exp){
    genelist <- read.csv('data/gene_signatures.csv',check.names = FALSE)
    genelist <- as.list(genelists)
    results <- gsva(exp,
                 genelist,
                 method = 'ssgsea',
                 ssgsea.norm=FALSE)
}

# run single sample gensetenrichment analysis on user defined geneset
ssGSEA_OV_custom <- function(exp,genelist_custom){
    genelist <- as.list(genelist_custom)
    results <- gsva(exp,
                 genelist,
                 method = 'ssgsea',
                 ssgsea.norm=FALSE)
}
