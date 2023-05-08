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
relassay <- assay(tse, "relabundance")

# Groundwork for plotting and K-means
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "bray")


# Simple hierarchical clustering on the rows
hclust.out <- clusterRows(relassay, HclustParam(), full=TRUE)
rowData(tse)$clusters <- hclust.out$clusters

# Dendogram plot for hierarchical clustering on rows (OMA book method)
dendro <- as.dendrogram(hclust.out$objects$hclust)
plot(dendro) # Classic dendogram
nbclusters <- length(levels(hclust.out$clusters))
col_val_map <- randomcoloR::distinctColorPalette(nbclusters) %>%
    as.list() %>% setNames(paste0("clust_",seq(nbclusters)))
dend <- color_branches(dendro, k=nbclusters, col=unlist(col_val_map))
labels(dend) <- NULL
plot(dend) # Plot with colors


# Simple hierarchical clustering on the cols
hclust.out <- clusterRows(t(assay(tse, "logcounts")), HclustParam(), full=TRUE)
colData(tse)$clusters <- hclust.out$clusters

# Plotting col clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")

# Dendogram plot for hierarchical clustering on cols (OMA book method)
tse <- runMDS(tse, assay.type = "logcounts", FUN = vegan::vegdist, method = "bray")

# Simple hierarchical clustering on the cols
dendro <- as.dendrogram(hclust.out$objects$hclust)
plot(dendro) # Classic dendogram
nbclusters <- length(levels(hclust.out$clusters))
col_val_map <- randomcoloR::distinctColorPalette(nbclusters) %>%
    as.list() %>% setNames(paste0("clust_",seq(nbclusters)))
dend <- color_branches(dendro, k=nbclusters, col=unlist(col_val_map))
labels(dend) <- NULL
plot(dend) # Plot with colors

# More complex hierarchical clustering on the rows: dynamic cut and different method
hclust.out <- clusterRows(assay(tse, "logcounts"), HclustParam(method="ward.D2", cut.dynamic=TRUE), full=TRUE)
hclust.out <- clusterRows(t(assay(tse, "logcounts")), HclustParam(method="ward.D2", cut.dynamic=TRUE), full=TRUE)
dendro <- as.dendrogram(hclust.out$objects$hclust)
plot(dendro) # Classic dendogram


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
kmeans.out <- clusterRows(relassay, KmeansParam(10))
head(kmeans.out, 10)
kmeans.out <- clusterRows(t(relassay), KmeansParam(10))
head(kmeans.out, 10)
plotReducedDim(tse, "MDS", colour_by = I(kmeans.out))


#DBSCAN
dbscan.out <- clusterRows(relassay, DbscanParam())
dbscan.out <- clusterRows(t(relassay), DbscanParam(eps=0.2, min.pts = 3))
colData(tse)$clusters <- dbscan.out
plotReducedDim(tse, "MDS", colour_by = "clusters")

# clusterSweep:
clusSweep <-clusterSweep(relassay, DbscanParam(), min.pts = c(1L, 2L, 3L, 4L), eps = c(0.05, 0.1, 0.15, 0.2))

# Self organizing maps
if( !require("kohonen") ){
    BiocManager::install("kohonen")
    library("kohonen")
}
set.seed(1000)

# On the rows
som.out <- clusterRows(relassay, SomParam(20))
head(som.out)

# On the cols
som.out <- clusterRows(t(relassay), SomParam(20), full=TRUE)
head(som.out$clusters)
colData(tse)$clusters <- som.out$clusters

#Plotting clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")

#Plotting relationship between clusters
par(mfrow=c(1,2))
plot(som.out$objects$som, "counts")
grid <- som.out$objects$som$grid$pts
text(grid[,1], grid[,2], seq_len(nrow(grid)))


# Graph based clustering
set.seed(101) # just in case there are ties.
graph.out <- clusterRows(relassay, NNGraphParam(k=10)) # On rows
head(graph.out, 10)
graph.out <- clusterRows(t(relassay), NNGraphParam(k=3), full=TRUE) # On cols
colData(tse)$clusters <- graph.out$clusters
head(graph.out, 10)
plotReducedDim(tse, "MDS", colour_by = "clusters")

# Two phase clustering 
set.seed(100) # for the k-means
two.out <- clusterRows(relassay, TwoStepParam(first = KmeansParam(centers=sqrt, iter.max=15))) # On rows
head(two.out, 50)

two.out <- clusterRows(t(relassay), TwoStepParam(first = KmeansParam(centers=sqrt, iter.max=10)), full=TRUE) # On cols
two.out$clusters
colData(tse)$clusters <- two.out$clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")


# k-medoids method
pam.out <- clusterRows(relassay, PamParam(centers=20), full = TRUE)
head(pam.out$clusters, 100)
pam.out <- clusterRows(t(relassay), PamParam(centers=4), full = TRUE)
colData(tse)$clusters <- pam.out$clusters
pam.out$clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")


# Mini-batch k-means method
if( !require("mbkmeans") ){
    BiocManager::install("mbkmeans")
    library("mbkmeans")
}
mbk.out <- clusterRows(relassay, bluster::MbkmeansParam(centers=20), full = TRUE)
head(mbk.out$clusters, 100)
mbk.out <- clusterRows(t(relassay), bluster::MbkmeansParam(centers=4), full = TRUE)
colData(tse)$clusters <- mbk.out$clusters
mbk.out$clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")

# Comparing algorithms
hclust.out <- clusterRows(relassay, HclustParam())
kmeans.out <- clusterRows(relassay, KmeansParam(10))
dbscan.out <- clusterRows(relassay, DbscanParam())
som.out <- clusterRows(relassay, SomParam(20))
graph.out <- clusterRows(relassay, NNGraphParam(k=10))
pam.out <- clusterRows(relassay, PamParam(centers=20))
mbk.out <- clusterRows(relassay, bluster::MbkmeansParam(centers=20))
compareClusterings(list(dbscan=dbscan.out, kmeans=kmeans.out, hclust=hclust.out, som=som.out, graph=graph.out, pam=pam.out, mbk=mbk.out))
# Result:
# dbscan        kmeans        hclust           som        graph          pam
# dbscan   0.000000000 -5.679752e-02 -2.568936e-03 -0.2807967607 1.467090e+00 1.301938e+02
# kmeans  -0.056797522  0.000000e+00  9.065185e-02  0.3894139727 3.315326e-05 9.369053e-02
# hclust  -0.002568936  9.065185e-02  0.000000e+00  0.0227030938 1.636843e-06 4.659129e-03
# som     -0.280796761  3.894140e-01  2.270309e-02  0.0000000000 1.560466e-04 3.323043e-01
# graph    1.467089745  3.315326e-05  1.636843e-06  0.0001560466 0.000000e+00 8.215992e-04
# pam    130.193757796  9.369053e-02  4.659129e-03  0.3323043048 8.215992e-04 0.000000e+00
# mbk     -1.689581187  1.507514e-01  7.711154e-03  0.3891957778 4.863923e-04 5.220366e-01
# mbk
# dbscan -1.6895811866
# kmeans  0.1507514328
# hclust  0.0077111545
# som     0.3891957778
# graph   0.0004863923
# pam     0.5220365797
# mbk     0.0000000000
