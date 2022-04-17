library(shiny)
library(ggdendro)
library(plotly)
library(dendextend)

# Setting up the env
library(WGCNA)
library(impute)
library(flashClust)
options(stringsAsFactors = FALSE)

enableWGCNAThreads()

WGCNAShinyUI <- function(id) {
	ns <- NS(id)
	
	tagList(h1("OutlierDetectionPlot"),
			plotlyOutput(ns("OutlierDetectionPlot")))
}

WGCNAShiny <- function(id) {
	moduleServer(id,
				 function(input, output, session) {
				 	# Load the data from the file
				 	# File downloaded from http://refex.dbcls.jp/download.php?lang=en
				 	path_to_file <-  here::here("RefEx_expression_EST10_human.tsv")
				 	datExpr <- read.delim(path_to_file, row.names = 1)
				 	# Replace -1 by NA
				 	datExpr[datExpr == '-1'] <- NA
				 	# View(datExpr)
				 	dim(datExpr)
				 	names(datExpr)
				 	
				 	
				 	datExpr = as.data.frame(lapply(datExpr, as.numeric))
				 	datExpr <- na.omit(datExpr)
				 	# Normalization with log2
				 	datExpr = log2(datExpr)
				 	
				 	head(datExpr[1:5, 1:5]) # samples in row, genes in column
				 	
				 	
				 	n = 1000
				 	
				 	datExpr0 = t(datExpr)[, 1:n]
				 	datExpr0 = datExpr[1:n, ]
				 	sampleTree = hclust(dist(datExpr0), method = "average")
				 	
				 	# sampleTree = hclust(dist(datExpr[1:n,]), method = "average");
				 	
				 	#####################
				 	dend2 <- sampleTree %>% as.dendrogram() %>% color_branches(3)
				 	p <- ggplot(dend2, mapping = aes_auto())
				 	# ggplotly(p)
				 	output$OutlierDetectionPlot <- renderPlotly(ggplotly(p))
				 	# output$OutlierDetectionPlot <- renderPlot(p)
				 	#####################
				 	
				 	cutHeight = 10
				 	# Plot a line to show the cut
				 	p <- p + geom_abline(slope = 0, intercept = cutHeight, col = "red")
				 	output$OutlierDetectionPlot <- renderPlotly(ggplotly(p))
				 	
				 	# Determine cluster under the line
				 	clust = cutreeStatic(sampleTree, cutHeight, minSize = ceiling(n * 0.05))
				 	table(clust)
				 	# clust 1 contains the samples we want to keep.
				 	keepSamples = (clust == 1)
				 	datExpr = datExpr0[keepSamples,]
				 	# traitColors = numbers2colors(datExpr, signed = FALSE, lim = c(min(datExpr, na.rm = TRUE), max(datExpr, na.rm = TRUE)));
				 	traitColors = numbers2colors(datExpr, signed = FALSE)
				 	
				 	
				 	# Re-cluster samples
				 	sampleTree2 = hclust(dist(datExpr), method = "average")
				 	# # Convert traits to a color representation: white means low, red means high, grey means missing entry
				 	# traitColors = numbers2colors(datExpr, signed = FALSE);
				 	# Plot the sample dendrogram and the colors underneath.
				 	plotDendroAndColors(
				 		sampleTree2,
				 		traitColors,
				 		groupLabels = names(datExpr),
				 		main = "Sample dendrogram and trait heatmap"
				 	)
				 	
				 	datExpr = t(datExpr)
				 	nGenes = ncol(datExpr)
				 	nSamples = nrow(datExpr)
				 	
				 	# Choose a set of soft-thresholding powers
				 	powers = c(c(1:10), seq(
				 		from = 12,
				 		to = 20,
				 		by = 2
				 	))
				 	# Call the network topology analysis function
				 	sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
				 	# sft = pickSoftThreshold(datExpr)
				 	# Plot the results:
				 	# sizeGrWindow(9, 5)
				 	par(mfrow = c(1, 2))
				 	
				 	cex1 = 0.9
				 	
				 	# Scale-free topology fit index as a function of the soft-thresholding power
				 	plot(
				 		sft$fitIndices[, 1],
				 		-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
				 		xlab = "Soft Threshold (power)",
				 		ylab = "Scale Free Topology Model Fit,signed R^2",
				 		type = "n",
				 		main = paste("Scale independence")
				 	)
				 	
				 	text(
				 		sft$fitIndices[, 1],
				 		-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
				 		labels = powers,
				 		cex = cex1,
				 		col = "red"
				 	)
				 	
				 	# this line corresponds to using an R^2 cut-off of h
				 	abline(h = 0.60, col = "red")
				 	# Mean connectivity as a function of the soft-thresholding power
				 	plot(
				 		sft$fitIndices[, 1],
				 		sft$fitIndices[, 5],
				 		xlab = "Soft Threshold (power)",
				 		ylab = "Mean Connectivity",
				 		type = "n",
				 		main = paste("Mean connectivity")
				 	)
				 	text(
				 		sft$fitIndices[, 1],
				 		sft$fitIndices[, 5],
				 		labels = powers,
				 		cex = cex1,
				 		col = "red"
				 	)
				 	net = blockwiseModules(
				 		datExpr,
				 		power = 10,
				 		deepSplit = 2,
				 		TOMType = "unsigned",
				 		minModuleSize = 30,
				 		reassignThreshold = 0,
				 		mergeCutHeight = 0.25,
				 		numericLabels = TRUE,
				 		pamRespectsDendro = FALSE,
				 		saveTOMs = TRUE,
				 		saveTOMFileBase = "humanTOM",
				 		verbose = 3
				 	)
				 	table(net$colors)
				 	# open a graphics window
				 	# sizeGrWindow(12, 9)
				 	# Convert labels to colors for plotting
				 	mergedColors = labels2colors(net$colors)
				 	# Plot the dendrogram and the module colors underneath
				 	plotDendroAndColors(
				 		net$dendrograms[[1]],
				 		mergedColors[net$blockGenes[[1]]],
				 		"Module colors",
				 		dendroLabels = FALSE,
				 		hang = 0.03,
				 		addGuide = TRUE,
				 		guideHang = 0.05
				 	)
				 	moduleLabels = net$colors
				 	moduleColors = labels2colors(net$colors)
				 	MEs = net$MEs
				 	
				 	geneTree = net$dendrograms[[1]]
				 	
				 	
				 	# Define numbers of genes and samples
				 	nGenes = ncol(datExpr)
				 	
				 	nSamples = nrow(datExpr)
				 	
				 	# Recalculate MEs with color labels
				 	MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
				 	MEs = orderMEs(MEs0)
				 	
				 	
				 	# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
				 	# calculated during module detection, but let us do it again here.
				 	dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = 10)
				 	
				 	# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
				 	plotTOM = dissTOM ^ 7
				 	
				 	# Set diagonal to NA for a nicer plot
				 	diag(plotTOM) = NA
				 	
				 	# Call the plot function
				 	# sizeGrWindow(9,9)
				 	TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
				 })
}
