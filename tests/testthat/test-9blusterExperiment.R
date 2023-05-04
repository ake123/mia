if( !require("scRNAseq") ){
    BiocManager::install("scRNAseq")
    library("scRNAseq")
}

if( !require("scran") ){
    BiocManager::install("scran")
    library("scran")
}


if( !require("scater") ){
    BiocManager::install("scater")
    library("scater")
}

#Data load
sce <- ZeiselBrainData()

# Trusting the authors' quality control, and going straight to normalization.
library(scuttle)
sce <- logNormCounts(sce)

# Feature selection based on highly variable genes.
library(scran)
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n=1000)

# Dimensionality reduction for work (PCA) and pleasure (t-SNE).
set.seed(1000)
library(scater)
sce <- runPCA(sce, ncomponents=20, subset_row=hvgs)
sce <- runUMAP(sce, dimred="PCA")

mat <- reducedDim(sce, "PCA")
dim(mat)

library(bluster)
hclust.out <- clusterRows(mat, HclustParam())
plotUMAP(sce, colour_by=I(hclust.out))
