# Package dependencies
if( !require("bluster") ){
    BiocManager::install("bluster")
    library("bluster")
}

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

if( !require("scuttle") ){
    BiocManager::install("scuttle")
    library("scuttle")
}


if( !require("dynamicTreeCut") ){
    BiocManager::install("dynamicTreeCut", INSTALL_opts = '--no-lock')
    library("dynamicTreeCut")
}

devtools::load_all()

# Data load
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns

# Normalization.
tse <- logNormCounts(tse)

# Clustering
hclust.out <- clusterRows(assay(tse, "logcounts"), HclustParam())
rowData(tse)$Cluster <- hclust.out



# Feature selection based on highly variable genes.
# dec <- modelGeneVar(tse)
# hvgs <- getTopHVGs(dec, n=1000)

# Dimensionality reduction for work (PCA) and pleasure (t-SNE).
set.seed(1000)
# tsetmp <- runPCA(tse, subset_row=rownames(tse), approx=FALSE)
tse <- runUMAP(tse, transposed=TRUE)

mat <- reducedDim(tsetmp, "PCA")
#dim(mat)

plotUMAP(tse, colour_by=I(hclust.out))

hclust.out <- clusterRows(assay(tse, "logcounts"), HclustParam(method="ward.D2", cut.dynamic=TRUE))


