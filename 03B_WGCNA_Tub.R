library(tidyverse)
library(WGCNA)
library(lmerTest)

####Load data####

load("./KI-DN_clin_hist.RData")

load("./KI_DN_clinBL_transcriptomics.RData")

####Data preparation of Tub data####

tub.meta <- dn.meta[dn.meta$Compartment=="Tub",]
tub.meta <- tub.meta[order(tub.meta$Exiqon.ID),]
tub.cases <- which(dn.meta$Compartment=="Tub")
tub.counts <- dn.counts[,c(tub.cases,79)]

tub.cases <- which(dn.meta$Status=="Tub-DN")
tub.counts.ckd <- dn.counts[,c(tub.cases,79)]
tub.meta.ckd <- dn.meta[dn.meta$Status=="Tub-DN",]

tub.meta <- merge(tub.meta, hist, by.x="Individual_ID", by.y="ID")
tub.meta <- merge(tub.meta, clinvar, by.x="Individual_ID", by.y="Pat.ID")
tub.meta.ckd <- tub.meta[tub.meta$Status=="Tub-DN",]
tub.meta.ckd <- tub.meta.ckd[order(tub.meta.ckd$Exiqon.ID),]

tub.counts.ckd.filt <- tub.counts.ckd
tub.counts.ckd.filt <- tub.counts.ckd.filt[apply(tub.counts.ckd.filt[,c(1:19)], 1, function(row) all(row !=0 )), ]

tub.tpm.ckd.filt <- dn.tpm[,c(tub.cases,79)]
tub.tpm.ckd.filt <- tub.tpm.ckd.filt[apply(tub.tpm.ckd.filt[,c(1:19)], 1, function(row) all(row >=1 )), ]

tpm.pass <- rownames(tub.tpm.ckd.filt)
tub.counts.ckd.filt <- tub.counts.ckd.filt[rownames(tub.counts.ckd.filt) %in% tpm.pass, ]

tpm.cts.pass <- rownames(tub.counts.ckd.filt)
tub.tpm.ckd.filt <- tub.tpm.ckd.filt[rownames(tub.tpm.ckd.filt) %in% tpm.cts.pass, ]

####Calculate eGFR slope####

df2 <- tub.meta.ckd %>% dplyr::select(Individual_ID, eGFR_FU1y, eGFR_FU2y, eGFR_FU3y, eGFR_FU4y, eGFR_FU5y)

df2.long <- df2 %>% 
  pivot_longer(!Individual_ID, names_to="Year", values_to = "eGFR")

df2.long <- df2.long %>% 
  mutate(Year = case_when(Year == "eGFR_FU1y" ~ 1,
                          Year == "eGFR_FU2y" ~ 2,
                          Year == "eGFR_FU3y" ~ 3,
                          Year == "eGFR_FU4y" ~ 4,
                          Year == "eGFR_FU5y" ~ 5))

baseline <- tub.meta.ckd %>% dplyr::select(Individual_ID, Gender, eGFR_BL, Diabetes_type)
baseline$Diabetes_type <- as.factor(baseline$Diabetes_type)

df2.long <- merge(df2.long, baseline, by.x="Individual_ID")


fm52 <- lmerTest::lmer(eGFR ~ Year + eGFR_BL + (0 + Year|Individual_ID), df2.long)
coef.fm52 <- as.data.frame(coef(fm52)$Individual_ID)

tub.meta.ckd <- merge(tub.meta.ckd, coef.fm52, by.x="Individual_ID", by.y="row.names")
tub.meta.ckd <- dplyr::rename(tub.meta.ckd, eGFR_slope = Year)
tub.meta.ckd$eGFR_BL.y <- NULL
tub.meta.ckd$`(Intercept)` <- NULL
tub.meta.ckd <- tub.meta.ckd[order(tub.meta.ckd$Exiqon.ID),]


####Run WGCNA####

wgcnaExpr = as.data.frame(t(tub.tpm.ckd.filt[, c(1:19)]))
names(wgcnaExpr) = tub.tpm.ckd.filt$symbol

sampleTree = fastcluster::hclust(dist(wgcnaExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

meta <- tub.meta.ckd[,c(1,3,9,11,12,21,83,84,140)]


meta$Exiqon.ID <- str_replace(meta$Exiqon.ID, "-", ".")
meta$Exiqon.ID <- str_c(meta$Exiqon.ID, "R", sep=".")
meta$Exiqon.ID <- str_c("X", meta$Exiqon.ID, sep="")

rownames(meta) = meta$Exiqon.ID
meta$Exiqon.ID <- NULL
meta$Individual_ID <- NULL
meta$Center <- NULL
meta$Gender <- NULL

wgcnaExpr <- wgcnaExpr[-c(2,6),]


# Re-cluster samples
sampleTree2 = fastcluster::hclust(dist(wgcnaExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(meta, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(meta), 
                    main = "Sample dendrogram and trait heatmap")


save(wgcnaExpr, meta, file = "WGCNATUB-01-dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=100, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(wgcnaExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net = blockwiseModules(wgcnaExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "WGCNATUB", 
                       verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
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
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "WGCNATUB-02-networkConstruction-auto.RData")


# Define numbers of genes and samples
nGenes = ncol(wgcnaExpr);
nSamples = nrow(wgcnaExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(wgcnaExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, meta, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(meta),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the weight column of datTrait
slope = as.data.frame(meta$eGFR_slope);
names(slope) = "eGFR slope"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(wgcnaExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(wgcnaExpr, slope, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(slope), sep="");
names(GSPvalue) = paste("p.GS.", names(slope), sep="");

# Create the starting data frame
geneInfo0 = data.frame(moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, slope, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.eGFR.slope));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "latest_geneInfotub.csv")



####Extract networks####


load(file = "WGCNATUB-01-dataInput.RData"); 
load(file = "WGCNATUB-02-networkConstruction-auto.RData");
softPower <- 9 ;
adjacency <- adjacency(wgcnaExpr, power = softPower) ;
TOM <- TOMsimilarity(adjacency) ;

tub.degs <- read_csv("LATEST_tub_degs.csv")

#BLACK

modules = c("black")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.black <- read_delim("CytoscapeInputtub-edges0-black.txt")
tub.edge.black %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.black
tub.edge.black %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.black

#yellow

modules = c("yellow")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.yellow <- read_delim("CytoscapeInputtub-edges0-yellow.txt")
tub.edge.yellow %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.yellow
tub.edge.yellow %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.yellow

#cyan

modules = c("cyan")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.cyan <- read_delim("CytoscapeInputtub-edges0-cyan.txt")
tub.edge.cyan %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.cyan
tub.edge.cyan %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.cyan

#royalblue

modules = c("royalblue")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.royalblue <- read_delim("CytoscapeInputtub-edges0-royalblue.txt")
tub.edge.royalblue %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.royalblue
tub.edge.royalblue %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.royalblue

#blue

modules = c("blue")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.blue <- read_delim("CytoscapeInputtub-edges0-blue.txt")
tub.edge.blue %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.blue
tub.edge.blue %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.blue


#greenyellow

modules = c("greenyellow")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.greenyellow <- read_delim("CytoscapeInputtub-edges0-greenyellow.txt")
tub.edge.greenyellow %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.greenyellow
tub.edge.greenyellow %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.greenyellow

#pink

modules = c("pink")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.pink <- read_delim("CytoscapeInputtub-edges0-pink.txt")
tub.edge.pink %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.pink
tub.edge.pink %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.pink

#grey60

modules = c("grey60")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.grey60 <- read_delim("CytoscapeInputtub-edges0-grey60.txt")
tub.edge.grey60 %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.grey60
tub.edge.grey60 %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.grey60

#salmon

modules = c("salmon")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.salmon <- read_delim("CytoscapeInputtub-edges0-salmon.txt")
tub.edge.salmon %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.salmon
tub.edge.salmon %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.salmon

#lightgreen

modules = c("lightgreen")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.lightgreen <- read_delim("CytoscapeInputtub-edges0-lightgreen.txt")
tub.edge.lightgreen %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.lightgreen
tub.edge.lightgreen %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.lightgreen

#green

modules = c("green")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.green <- read_delim("CytoscapeInputtub-edges0-green.txt")
tub.edge.green %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.green
tub.edge.green %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.green

#magenta

modules = c("magenta")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.magenta <- read_delim("CytoscapeInputtub-edges0-magenta.txt")
tub.edge.magenta %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.magenta
tub.edge.magenta %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.magenta

#turquoise

modules = c("turquoise")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.turquoise <- read_delim("CytoscapeInputtub-edges0-turquoise.txt")
tub.edge.turquoise %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.turquoise
tub.edge.turquoise %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.turquoise

#lightcyan

modules = c("lightcyan")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.lightcyan <- read_delim("CytoscapeInputtub-edges0-lightcyan.txt")
tub.edge.lightcyan %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.lightcyan
tub.edge.lightcyan %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.lightcyan

#midnightblue

modules = c("midnightblue")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.midnightblue <- read_delim("CytoscapeInputtub-edges0-midnightblue.txt")
tub.edge.midnightblue %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.midnightblue
tub.edge.midnightblue %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.midnightblue

#red

modules = c("red")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.red <- read_delim("CytoscapeInputtub-edges0-red.txt")
tub.edge.red %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.red
tub.edge.red %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.red

#lightyellow

modules = c("lightyellow")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.lightyellow <- read_delim("CytoscapeInputtub-edges0-lightyellow.txt")
tub.edge.lightyellow %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.lightyellow
tub.edge.lightyellow %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.lightyellow

#brown

modules = c("brown")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.brown <- read_delim("CytoscapeInputtub-edges0-brown.txt")
tub.edge.brown %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.brown
tub.edge.brown %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.brown

#purple

modules = c("purple")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.purple <- read_delim("CytoscapeInputtub-edges0-purple.txt")
tub.edge.purple %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.purple
tub.edge.purple %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.purple

#tan

modules = c("tan")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.tan <- read_delim("CytoscapeInputtub-edges0-tan.txt")
tub.edge.tan %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.tan
tub.edge.tan %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.tan

#grey

modules = c("grey")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputtub-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputtub-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

tub.edge.grey <- read_delim("CytoscapeInputtub-edges0-grey.txt")
tub.edge.grey %>% dplyr::filter(fromNode %in% tub.degs$symbol) -> tub.edge.grey
tub.edge.grey %>% dplyr::filter(toNode %in% tub.degs$symbol) -> tub.edge.grey


rbind(tub.edge.black, tub.edge.yellow, tub.edge.cyan, tub.edge.royalblue,
      tub.edge.blue, tub.edge.pink, tub.edge.salmon, tub.edge.lightgreen,
      tub.edge.green, tub.edge.turquoise, tub.edge.red, tub.edge.grey) -> wgcna.tub.network

writexl::write_xlsx(wgcna.tub.network, "wgcna_tub_network.xlsx")

