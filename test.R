# Setting up the env
library(WGCNA)
library(impute)
library(flashClust)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Load the data from the file
# File downloaded from http://refex.dbcls.jp/download.php?lang=en
path_to_file <-  here::here("RefEx_expression_EST10_human.tsv")
datExpr <- read.delim(path_to_file, row.names = 1)
# Replace -1 by NA
datExpr[datExpr == '-1'] <- NA

# View(datExpr)
dim(datExpr)
names(datExpr)
# datExpr0 = as.data.frame(t(datExpr[, -c(1:8)]));
# names(datExpr0) = datExpr$substanceBXH;
# rownames(datExpr0) = names(datExpr)[-c(1:8)];
# View(datExpr0)
# View(datExpr)
# gsg = goodSamplesGenes(datExpr0, verbose = 3);
# gsg$allOK

datExpr = as.data.frame(lapply(datExpr, as.numeric))
datExpr <- na.omit(datExpr)
# Normalization with log2
datExpr = log2(datExpr)

head(datExpr[1:5,1:5]) # samples in row, genes in column

n = 5000;
datExpr0=t(datExpr)[,1:n]
datExpr0=datExpr[1:n,]
sampleTree = hclust(dist(datExpr0), method = "average");
# sampleTree = hclust(dist(datExpr[1:n,]), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
	 cex.axis = 1.5, cex.main = 2)
cutHeight = 10
# Plot a line to show the cut
abline(h = cutHeight, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight, minSize = ceiling(n*0.05))
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

# traitColors = numbers2colors(datExpr, signed = FALSE, lim = c(min(datExpr, na.rm = TRUE), max(datExpr, na.rm = TRUE)));
traitColors = numbers2colors(datExpr, signed = FALSE);
# traitData = read.csv("ClinicalTraits.csv");
# dim(traitData)
# names(traitData)
# # remove columns that hold information we do not need.
# allTraits = traitData[, -c(31, 16)];
# allTraits = allTraits[, c(2, 11:36) ];
# dim(allTraits)
# names(allTraits)
# # Form a data frame analogous to expression data that will hold the clinical traits.
# femaleSamples = rownames(datExpr);
# traitRows = match(femaleSamples, allTraits$Mice);
# datTraits = allTraits[traitRows, -1];
# rownames(datTraits) = allTraits[traitRows, 1];
# collectGarbage();
# # Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# # Convert traits to a color representation: white means low, red means high, grey means missing entry
# traitColors = numbers2colors(datExpr, signed = FALSE);
# # Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
					groupLabels = names(datExpr),
					main = "Sample dendrogram and trait heatmap")

datExpr = t(datExpr)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# sft = pickSoftThreshold(datExpr)
# Plot the results:
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	 xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	 main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	 labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.50,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
	 xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	 main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
net = blockwiseModules(datExpr, power = 10, deepSplit=2,
					   TOMType = "unsigned", minModuleSize = 30,
					   reassignThreshold = 0, mergeCutHeight = 0.25,
					   numericLabels = TRUE, pamRespectsDendro = FALSE,
					   saveTOMs = TRUE,
					   saveTOMFileBase = "humanTOM",
					   verbose = 3)
table(net$colors)
# open a graphics window
# sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
					"Module colors",
					dendroLabels = FALSE, hang = 0.03,
					addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# moduleTraitCor = cor(MEs, datExpr, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# sizeGrWindow(10,6)
# # Will display correlations and their p-values
# textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
# signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,
# 			   xLabels = names(t(datExpr)),
# 			   yLabels = names(MEs),
# 			   ySymbols = names(MEs),
# 			   colorLabels = FALSE,
# 			   colors = greenWhiteRed(50),
# 			   textMatrix = textMatrix,
# 			   setStdMargins = FALSE,
# 			   cex.text = 0.5,
# 			   zlim = c(-1,1),
# 			   main = paste("Module-trait relationships"))
# # sizeGrWindow(10,6)
# # Will display correlations and their p-values
# textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
# 				   signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,
# 			   xLabels = names(datTraits),
# 			   yLabels = names(MEs),
# 			   ySymbols = names(MEs),
# 			   colorLabels = FALSE,
# 			   colors = blueWhiteRed(50),
# 			   textMatrix = textMatrix,
# 			   setStdMargins = FALSE,
# 			   cex.text = 0.5,
# 			   zlim = c(-1,1),
# 			   main = paste("Module-trait relationships"))
# # Define variable weight containing the weight column of datTrait
# weight = as.data.frame(datTraits$weight_g);
# names(weight) = "weight"
# # names (colors) of the modules
# modNames = substring(names(MEs), 3)
# geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
# names(geneModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p.MM", modNames, sep="");
# geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
# names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
# names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# module = "brown"
# column = match(module, modNames);
# moduleGenes = moduleColors==module;
# # sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
# 				   abs(geneTraitSignificance[moduleGenes, 1]),
# 				   xlab = paste("Module Membership in", module, "module"),
# 				   ylab = "Gene significance for body weight",
# 				   main = paste("Module membership vs. gene significance\n"),
# 				   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# names(datExpr)
# annot = read.csv(file = "GeneAnnotation.csv");
# dim(annot)
# names(annot)
# probes = names(datExpr)
# probes2annot = match(probes, annot$substanceBXH)
# # The following is the number or probes without annotation:
# sum(is.na(probes2annot))
# # Should return 0.
# # Create the starting data frame
# geneInfo0 = data.frame(substanceBXH = probes,
# 					   geneSymbol = annot$gene_symbol[probes2annot],
# 					   LocusLinkID = annot$LocusLinkID[probes2annot],
# 					   moduleColor = moduleColors,
# 					   geneTraitSignificance,
# 					   GSPvalue)
# # Order modules by their significance for weight
# modOrder = order(-abs(cor(MEs, weight, use = "p")));
# # Add module membership information in the chosen order
# for (mod in 1:ncol(geneModuleMembership))
# {
# 	oldNames = names(geneInfo0)
# 	geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
# 						   MMPvalue[, modOrder[mod]]);
# 	names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
# 						 paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }
# # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
# geneInfo = geneInfo0[geneOrder, ]
# # Create the starting data frame
# geneInfo0 = data.frame(substanceBXH = probes,
# 					   geneSymbol = annot$gene_symbol[probes2annot],
# 					   LocusLinkID = annot$LocusLinkID[probes2annot],
# 					   moduleColor = moduleColors,
# 					   geneTraitSignificance,
# 					   GSPvalue)
# # Order modules by their significance for weight
# modOrder = order(-abs(cor(MEs, weight, use = "p")));
# # Add module membership information in the chosen order
# for (mod in 1:ncol(geneModuleMembership))
# {
# 	oldNames = names(geneInfo0)
# 	geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
# 						   MMPvalue[, modOrder[mod]]);
# 	names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
# 						 paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }
# # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
# geneInfo = geneInfo0[geneOrder, ]
# write.csv(geneInfo, file = "geneInfo.csv")
# fix(geneInfo)
# View(geneInfo)

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 10);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
# sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

# # Recalculate module eigengenes
# MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# # Isolate weight from the clinical traits
# weight = as.data.frame(datTraits$weight_g);
# names(weight) = "weight"
# # Add the weight to existing module eigengenes
# MET = orderMEs(cbind(MEs, weight))
# # Plot the relationships among the eigengenes and the trait
# # sizeGrWindow(5,7.5);
# par(cex = 0.9)
# plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
# 					  = 90)
# # Plot the dendrogram
# # sizeGrWindow(6,6);
# par(cex = 1.0)
# plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
# 					  plotHeatmaps = FALSE)
# # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
# par(cex = 1.0)
# plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
# 					  plotDendrograms = FALSE, xLabelsAngle = 90)
# dist(datExpr)
# # Plot the dendrogram
# # sizeGrWindow(6,6);
# par(cex = 1.0)
# plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
# 					  plotHeatmaps = FALSE)
# MET = orderMEs(cbind(MEs, weight))
# # Plot the relationships among the eigengenes and the trait
# # sizeGrWindow(5,7.5);
# par(cex = 0.9)
# plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
# 					  = 90)

# install.packages("reprex")
# install.packages(c("commonmark", "desc", "Rcpp", "rlang", "rmarkdown", "survival"))
# savehistory("D:/Programs/WGCNA_Trial/temp.r")
