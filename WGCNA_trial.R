# Setting up the env
library(WGCNA)
library(impute)
library(flashClust)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Load the data from the file
# path_to_file <-  here::here("RefEx_expression_EST10_human.tsv")
# Database <- read.delim(path_to_file, row.names = 1)

# load("oed.RData")
# Database <- oed

dim(Database)

# Following steps from http://pklab.med.harvard.edu/scw2014/WGCNA.html
gene.names=rownames(Database)

#Take a subset of genes, WGCNA requires genes be given in the columns
n=500;
datExpr=t(Database)[,1:n]
dim(datExpr)

SubGeneNames=gene.names[1:n]


# Remove NAs
datExpr <- na.omit(datExpr)
dim(datExpr)

# Choosing a soft-threshold to fit a scale-free topology to the network
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

# Plots the result
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# Generating adjacency and TOM similarity matrices based on the selected softpower
softPower = 7

# Calculates the adjacency matrix
adj = adjacency(datExpr,type = "unsigned", power = softPower)

# Turns adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM = TOMsimilarityFromExpr(datExpr)

colnames(TOM) = rownames(TOM) = SubGeneNames
dissTOM = 1-TOM


# Module detection

# hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM), method="average")

# plots the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", cex=0.3)

