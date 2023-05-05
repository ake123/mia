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

if( !require("apcluster") ){
    BiocManager::install("apcluster", INSTALL_opts = '--no-lock')
    library("apcluster")
}

if( !require(dendextend) ){
    install.packages("dendextend")
    library(dendextend)
}

devtools::load_all()

# Data load
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns

# Normalization.
tse <- logNormCounts(tse)
tse <- transformCounts(tse, method = "relabundance")

# Groundwork for plotting
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "bray")
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "bray", transposed=T)

# Simple hierarchical clustering on the rows
hclust.out <- clusterRows(assay(tse, "logcounts"), HclustParam(), full=TRUE)
rowData(tse)$clusters <- hclust.out$clusters

# Simple hierarchical clustering on the cols
hclust.out <- clusterRows(t(assay(tse, "logcounts")), HclustParam(), full=TRUE)
colData(tse)$clusters <- hclust.out$clusters

# Plotting clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")

#Plotting dendogram with OMA Book method
dendro <- as.dendrogram(hclust.out$objects$hclust)
plot(dendro) # Classic dendogram
nbclusters <- length(levels(hclust.out$clusters))
col_val_map <- randomcoloR::distinctColorPalette(nbclusters) %>%
    as.list() %>% setNames(paste0("clust_",seq(nbclusters)))
dend <- color_branches(dendro, k=nbclusters, col=unlist(col_val_map))
labels(dend) <- NULL
plot(dend) # Plot with colors
tse <- runMDS(tse, assay.type = "logcounts", FUN = vegan::vegdist, method = "bray")

# Simple hierarchical clustering on the cols
hclust.out <- clusterRows(t(assay(tse, "logcounts")), HclustParam(), full=TRUE)
colData(tse)$clusters <- hclust.out$clusters
dendro <- as.dendrogram(hclust.out$objects$hclust)
plot(dendro) # Classic dendogram
nbclusters <- length(levels(hclust.out$clusters))
col_val_map <- randomcoloR::distinctColorPalette(nbclusters) %>%
    as.list() %>% setNames(paste0("clust_",seq(nbclusters)))
dend <- color_branches(dendro, k=nbclusters, col=unlist(col_val_map))
labels(dend) <- NULL
plot(dend) # Plot with colors
# More complex hierarchical clustering on the rows: dynamic cut and different method
hclust.out <- clusterRows(assay(tse, "logcounts"), HclustParam(method="ward.D2", cut.dynamic=TRUE))
hclust.out <- clusterRows(t(assay(tse, "logcounts")), HclustParam(method="ward.D2", cut.dynamic=TRUE))


# Affinity propagation: quite slow
set.seed(1000)
sub <- tse[sample(nrow(tse), 1000),sample(ncol(tse), 26)]
# By row
ap.out <- clusterRows(assay(sub, "logcounts"), AffinityParam(), full=T)

#By column
ap.out <- clusterRows(t(assay(sub, "logcounts")), AffinityParam(), full=T)
colData(tse)$clusters <- ap.out$clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")


# K-means clustering
set.seed(100)

# By row
kmeans.out <- clusterRows(reducedDim(tse), KmeansParam(10))
plotReducedDim(tse, "MDS", colour_by = I(kmeans.out))
