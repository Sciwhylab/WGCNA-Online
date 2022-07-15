library(shiny)
# library(ggdendro)
# library(dendextend)

# Setting up the env
library(WGCNA)
library(impute)
library(flashClust)
options(stringsAsFactors = FALSE)

enableWGCNAThreads()

WGCNAShinyUI <- function(id) {
	ns <- NS(id)
	
	tagList(
		h1("Analysis"),
		
		sliderInput(
			ns("NoOfSamples"),
			label = "Sample Size",
			value = 1000,
			min = 100,
			max = 10000,
			step = 100,
			width = "100%"
		),
		
		h2("Outlier Detection"),
		numericInput(ns("cutHeight"), label = "Cut Height", value = 10),
		plotOutput(ns("OutlierDetectionPlot")),
		h2("Sample dendrogram and trait heatmap"),
		plotOutput(ns("DendroAndColorsPlot")),
		h2("Soft-thresholding"),
		# numericInput(
		# 	ns("h"),
		# 	label = "h for R^2 cut-off",
		# 	value = 0.1,
		# 	step = 0.05
		# ),
		selectInput(
			ns("power"),
			label = "Power for R^2 cut-off",
			choices = seq.int(from = 1, to = 30, by = 1)
		),
		plotOutput(ns("ScaleIndependencePlot")),
		plotOutput(ns("MeanConnectivityPlot")),
		h2("Dendrogram and the Module Colours"),
		plotOutput(ns("ColouredDendrogramPlot")),
		h2("Network heatmap plot, all genes"),
		plotOutput(ns("NetworkHeatmapPlot"))
	)
}

WGCNAShiny <- function(id, Dataset) {
	moduleServer(id,
				 function(input, output, session, datExpr1=Dataset) {
				 	
				 	n = reactive({
				 		input$NoOfSamples
				 	})
				 	
				 	# datExpr0 <- reactive({
				 	# 	t(datExpr1)[, 1:n()]
				 	# })
				 	datExpr0 <- reactive({
				 		datExpr1()[1:n(), ]
				 	})
				 	sampleTree <-
				 		reactive({
				 			hclust(dist(datExpr0()), method = "average")
				 		})
				 	
				 	# sampleTree = hclust(dist(datExpr[1:n(),]), method = "average");
				 	
				 	#####################
				 	# dend <- sampleTree %>% as.dendrogram() %>% color_branches(3)
				 	# ggplotly(p)
				 	# output$OutlierDetectionPlot <- renderPlot(p)
				 	#####################
				 	
				 	cutHeight <- reactive({
				 		input$cutHeight
				 	})
				 	output$OutlierDetectionPlot <- renderCachedPlot({
				 		par(cex = 0.6)
				 		
				 		par(mar = c(0, 4, 2, 0))
				 		plot(
				 			sampleTree(),
				 			main = "Sample clustering to detect outliers",
				 			sub = "",
				 			xlab = "",
				 			cex.lab = 1.5,
				 			cex.axis = 1.5,
				 			cex.main = 2
				 		)
				 		# Plot a line to show the cut
				 		abline(h = cutHeight(), col = "red")
				 	},
				 	cacheKeyExpr = {
				 		list(sampleTree(), cutHeight())
				 	})
				 	
				 	# Determine cluster under the line
				 	clust <-
				 		reactive({
				 			cutreeStatic(sampleTree(), cutHeight(), minSize = ceiling(n() * 0.05))
				 		})
				 	# table(clust())
				 	# clust 1 contains the samples we want to keep.
				 	keepSamples <- reactive({
				 		(clust() == 1)
				 	})
				 	datExpr <- reactive({
				 		datExpr0()[keepSamples(), ]
				 	})
				 	# traitColors = numbers2colors(datExpr, signed = FALSE, lim = c(min(datExpr, na.rm = TRUE), max(datExpr, na.rm = TRUE)));
				 	# Convert traits to a color representation: white means low, red means high, grey means missing entry
				 	traitColors = reactive({
				 		numbers2colors(datExpr(), signed = FALSE)
				 	})
				 	
				 	
				 	# Re-cluster samples
				 	sampleTree2 = reactive({
				 		hclust(dist(datExpr()), method = "average")
				 	})
				 	# Plot the sample dendrogram and the colors underneath.
				 	output$DendroAndColorsPlot <- renderCachedPlot({
				 		plotDendroAndColors(
				 			sampleTree2(),
				 			traitColors(),
				 			groupLabels = names(datExpr()),
				 			main = "Sample dendrogram and trait heatmap"
				 		)
				 	},
				 	cacheKeyExpr = {
				 		list(sampleTree2(), traitColors(), names(datExpr()))
				 	})
				 	
				 	datExpr3 <- reactive({
				 		t(datExpr())
				 	})
				 	# datExpr <- reactive({
				 	# 	t(datExpr0()[keepSamples(), ])
				 	# })
				 	# nGenes = ncol(datExpr3)
				 	# nSamples = nrow(datExpr3)
				 	
				 	# Choose a set of soft-thresholding powers
				 	powers <- c(c(1:10), seq(
				 		from = 12,
				 		to = 20,
				 		by = 2
				 	))
				 	# Call the network topology analysis function
				 	sft = reactive({
				 		pickSoftThreshold(datExpr3(),
				 						  powerVector = powers,
				 						  verbose = 5)
				 	}) %>% bindCache(datExpr3(), powers)
				 	# sft = pickSoftThreshold(datExpr3)
				 	
				 	observe({
				 		updateSelectInput(session, "power",
				 						  choices = round(sft()$fitIndices[, 1]), )
				 		# choices = seq.int(from = 1, to = 30, by = 1))
				 	})
				 	
				 	# bindEvent({ # Debugging only
				 	# 	output$Debug <- renderText({
				 	# 		paste(power(), typeof(power()))
				 	# 	})
				 	# },
				 	# input$power)
				 	
				 	# Scale-free topology fit index as a function of the soft-thresholding power
				 	output$ScaleIndependencePlot <- renderCachedPlot({
				 		cex1 <- 0.9
				 		par(mfrow = c(1, 2))
				 		
				 		plot(
				 			sft()$fitIndices[, 1],-sign(sft()$fitIndices[, 3]) * sft()$fitIndices[, 2],
				 			xlab = "Soft Threshold (power)",
				 			ylab = "Scale Free Topology Model Fit,signed R^2",
				 			type = "n",
				 			main = paste("Scale independence")
				 		)
				 		
				 		text(
				 			sft()$fitIndices[, 1],-sign(sft()$fitIndices[, 3]) * sft()$fitIndices[, 2],
				 			labels = powers,
				 			cex = cex1,
				 			col = "red"
				 		)
				 		
				 		# this line corresponds to using an R^2 cut-off of h
				 		abline(h = (-sign(sft()$fitIndices[, 3]) * sft()$fitIndices[, 2])[which(sft()$fitIndices[, 1] == power())[1]],
				 			   col = "red")
				 		
				 	},
				 	cacheKeyExpr = {
				 		list(sft(), power(), powers, input$h)
				 	})
				 	
				 	# Mean connectivity as a function of the soft-thresholding power
				 	output$MeanConnectivityPlot <- renderCachedPlot({
				 		cex1 <- 0.9
				 		par(mfrow = c(1, 2))
				 		
				 		plot(
				 			sft()$fitIndices[, 1],
				 			sft()$fitIndices[, 5],
				 			xlab = "Soft Threshold (power)",
				 			ylab = "Mean Connectivity",
				 			type = "n",
				 			main = paste("Mean connectivity")
				 		)
				 		text(
				 			sft()$fitIndices[, 1],
				 			sft()$fitIndices[, 5],
				 			labels = powers,
				 			cex = cex1,
				 			col = "red"
				 		)
				 	},
				 	cacheKeyExpr = {
				 		list(sft(), powers)
				 	})
				 	
				 	power <- reactive({
				 		as.double(input$power)
				 	})
				 	net = reactive({
				 		blockwiseModules(
				 			datExpr3(),
				 			power = power(),
				 			deepSplit = 2,
				 			TOMType = "unsigned",
				 			minModuleSize = 10,
				 			reassignThreshold = 0,
				 			mergeCutHeight = 0.25,
				 			numericLabels = TRUE,
				 			pamRespectsDendro = FALSE,
				 			# saveTOMs = TRUE,
				 			# saveTOMFileBase = "humanTOM",
				 			verbose = 3
				 		)
				 	}) %>%
				 		bindCache(datExpr3(), power())
				 	
				 	# table(net$colors)
				 	# # open a graphics window
				 	# # sizeGrWindow(12, 9)
				 	# Convert labels to colors for plotting
				 	mergedColors <- reactive({
				 		labels2colors(net()$colors)
				 	})
				 	# Plot the dendrogram and the module colours underneath
				 	output$ColouredDendrogramPlot <- renderCachedPlot({
				 		plotDendroAndColors(
				 			net()$dendrograms[[1]],
				 			mergedColors()[net()$blockGenes[[1]]],
				 			"Module colors",
				 			dendroLabels = FALSE,
				 			hang = 0.03,
				 			addGuide = TRUE,
				 			guideHang = 0.05
				 		)
				 	},
				 	cacheKeyExpr = {
				 		list(net(), mergedColors())
				 	})
				 	
				 	moduleLabels <- reactive({
				 		net()$colors
				 	})
				 	moduleColors <- reactive({
				 		labels2colors(net()$colors)
				 	})
				 	MEs <- reactive({
				 		net()$MEs
				 	})
				 	
				 	geneTree <- reactive({
				 		net()$dendrograms[[1]]
				 	})
				 	
				 	
				 	# # Define numbers of genes and samples
				 	# nGenes = ncol(datExpr3)
				 	#
				 	# nSamples = nrow(datExpr3)
				 	#
				 	# # Recalculate MEs with color labels
				 	# MEs0 = moduleEigengenes(datExpr3, moduleColors)$eigengenes
				 	# MEs = orderMEs(MEs0)
				 	
				 	
				 	# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
				 	# calculated during module detection, but let us do it again here.
				 	dissTOM <- reactive({
				 		1 - TOMsimilarityFromExpr(datExpr3(), power = power())
				 	}) %>%
				 		bindCache(datExpr3(), power())
				 	
				 	
				 	# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
				 	plotTOM <- reactive({
				 		dissTOM() ^ 7
				 	})
				 	
				 	# # Set diagonal to NA for a nicer plot
				 	# diag(plotTOM) = NA
				 	
				 	# Call the plot function
				 	output$NetworkHeatmapPlot <- renderCachedPlot({
				 		TOMplot(plotTOM(), geneTree(), moduleColors(), main = "Network heatmap plot, all genes")
				 	},
				 	cacheKeyExpr = {
				 		list(plotTOM(), geneTree(), moduleColors())
				 	})
				 	return(datExpr1)
				 })
}
