library(tidyverse)
library(WGCNA)
library(lmerTest)

####Load data####

load("./KI-DN_clin_hist.RData")

load("./KI_DN_clinBL_transcriptomics.RData")

####Data preparation of Glom data####

glom.meta <- dn.meta[dn.meta$Compartment=="Glom",]
glom.meta <- glom.meta[order(glom.meta$Exiqon.ID),]
glom.cases <- which(dn.meta$Compartment=="Glom")
glom.counts <- dn.counts[,c(glom.cases,79)]

glom.cases <- which(dn.meta$Status=="Glom-DN")
glom.counts.ckd <- dn.counts[,c(glom.cases,79)]
glom.meta.ckd <- dn.meta[dn.meta$Status=="Glom-DN",]

glom.meta <- merge(glom.meta, hist, by.x="Individual_ID", by.y="ID")
glom.meta <- merge(glom.meta, clinvar, by.x="Individual_ID", by.y="Pat.ID")
glom.meta.ckd <- glom.meta[glom.meta$Status=="Glom-DN",]
glom.meta.ckd <- glom.meta.ckd[order(glom.meta.ckd$Exiqon.ID),]

glom.counts.ckd.filt <- glom.counts.ckd
glom.counts.ckd.filt <- glom.counts.ckd.filt[apply(glom.counts.ckd.filt[,c(1:19)], 1, function(row) all(row !=0 )), ]

glom.tpm.ckd.filt <- dn.tpm[,c(glom.cases,79)]
glom.tpm.ckd.filt <- glom.tpm.ckd.filt[apply(glom.tpm.ckd.filt[,c(1:19)], 1, function(row) all(row >=1 )), ]

tpm.pass <- rownames(glom.tpm.ckd.filt)
glom.counts.ckd.filt <- glom.counts.ckd.filt[rownames(glom.counts.ckd.filt) %in% tpm.pass, ]

tpm.cts.pass <- rownames(glom.counts.ckd.filt)
glom.tpm.ckd.filt <- glom.tpm.ckd.filt[rownames(glom.tpm.ckd.filt) %in% tpm.cts.pass, ]

####Calculate eGFR slope####

df2 <- glom.meta.ckd %>% dplyr::select(Individual_ID, eGFR_FU1y, eGFR_FU2y, eGFR_FU3y, eGFR_FU4y, eGFR_FU5y)

df2.long <- df2 %>% 
  pivot_longer(!Individual_ID, names_to="Year", values_to = "eGFR")

df2.long <- df2.long %>% 
  mutate(Year = case_when(Year == "eGFR_FU1y" ~ 1,
                          Year == "eGFR_FU2y" ~ 2,
                          Year == "eGFR_FU3y" ~ 3,
                          Year == "eGFR_FU4y" ~ 4,
                          Year == "eGFR_FU5y" ~ 5))

baseline <- glom.meta.ckd %>% dplyr::select(Individual_ID, Gender, eGFR_BL, Diabetes_type)
baseline$Diabetes_type <- as.factor(baseline$Diabetes_type)

df2.long <- merge(df2.long, baseline, by.x="Individual_ID")


fm52 <- lmerTest::lmer(eGFR ~ Year + eGFR_BL + (0 + Year|Individual_ID), df2.long)
coef.fm52 <- as.data.frame(coef(fm52)$Individual_ID)

glom.meta.ckd <- merge(glom.meta.ckd, coef.fm52, by.x="Individual_ID", by.y="row.names")
glom.meta.ckd <- dplyr::rename(glom.meta.ckd, eGFR_slope = Year)
glom.meta.ckd$eGFR_BL.y <- NULL
glom.meta.ckd$`(Intercept)` <- NULL
glom.meta.ckd <- glom.meta.ckd[order(glom.meta.ckd$Exiqon.ID),]

glom.meta.ckd[7, 84] = 104.02

####Run WGCNA####

wgcnaExpr = as.data.frame(t(glom.tpm.batch.adj.df[, c(1:19)]))
names(wgcnaExpr) = glom.tpm.batch.adj.df$symbol

sampleTree = fastcluster::hclust(dist(wgcnaExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


meta <- glom.meta.ckd[,c(1,3,9,11,12,21,83,84,140)]


meta$Exiqon.ID <- str_replace(meta$Exiqon.ID, "-", ".")
meta$Exiqon.ID <- str_c(meta$Exiqon.ID, "R", sep=".")
meta$Exiqon.ID <- str_c("X", meta$Exiqon.ID, sep="")

rownames(meta) = meta$Exiqon.ID
meta$Exiqon.ID <- NULL
meta$Individual_ID <- NULL
meta$Center <- NULL
meta$Gender <- NULL

wgcnaExpr <- wgcnaExpr[-c(13,14),]


# Re-cluster samples
sampleTree2 = fastcluster::hclust(dist(wgcnaExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(meta, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(meta), 
                    main = "Sample dendrogram and trait heatmap")


save(wgcnaExpr, meta, file = "WGCNAGLOMadj-01-dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
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
                       saveTOMFileBase = "WGCNAGLOMadj", 
                       verbose = 3)


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
     file = "WGCNAGLOMadj-02-networkConstruction-auto.RData")


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



names(wgcnaExpr)

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


write.csv(geneInfo, file = "latest_geneInfoglomadj.csv")


####Extract networks####

load(file = "WGCNAGLOMadj-01-dataInput.RData"); 
load(file = "WGCNAGLOMadj-02-networkConstruction-auto.RData");
softPower <- 9 ;
adjacency <- adjacency(wgcnaExpr, power = softPower) ;
TOM <- TOMsimilarity(adjacency) ;


glom.degs <- read_csv("LATEST_glom_degs.csv")

#BLACK

modules = c("black")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.black <- read_delim("CytoscapeInputglom-edges0-black.txt")
glom.edge.black %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.black
glom.edge.black %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.black

#blue

modules = c("blue")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.blue <- read_delim("CytoscapeInputglom-edges0-blue.txt")
glom.edge.blue %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.blue
glom.edge.blue %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.blue

#brown

modules = c("brown")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.brown <- read_delim("CytoscapeInputglom-edges0-brown.txt")
glom.edge.brown %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.brown
glom.edge.brown %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.brown

#cyan

modules = c("cyan")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.cyan <- read_delim("CytoscapeInputglom-edges0-cyan.txt")
glom.edge.cyan %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.cyan
glom.edge.cyan %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.cyan

#darkgreen

modules = c("darkgreen")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);


glom.edge.darkgreen <- read_delim("CytoscapeInputglom-edges0-darkgreen.txt")
glom.edge.darkgreen %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.darkgreen
glom.edge.darkgreen %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.darkgreen

#darkgrey

modules = c("darkgrey")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.darkgrey <- read_delim("CytoscapeInputglom-edges0-darkgrey.txt")
glom.edge.darkgrey %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.darkgrey
glom.edge.darkgrey %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.darkgrey

#darkorange

modules = c("darkorange")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.darkorange <- read_delim("CytoscapeInputglom-edges0-darkorange.txt")
glom.edge.darkorange %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.darkorange
glom.edge.darkorange %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.darkorange

#darkred

modules = c("darkred")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);


glom.edge.darkred <- read_delim("CytoscapeInputglom-edges0-darkred.txt")
glom.edge.darkred %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.darkred
glom.edge.darkred %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.darkred

#darkturquoise

modules = c("darkturquoise")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);


glom.edge.darkturquoise <- read_delim("CytoscapeInputglom-edges0-darkturquoise.txt")
glom.edge.darkturquoise %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.darkturquoise
glom.edge.darkturquoise %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.darkturquoise

#green

modules = c("green")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.green <- read_delim("CytoscapeInputglom-edges0-green.txt")
glom.edge.green %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.green
glom.edge.green %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.green

#greenyellow

modules = c("greenyellow")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.greenyellow <- read_delim("CytoscapeInputglom-edges0-greenyellow.txt")
glom.edge.greenyellow %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.greenyellow
glom.edge.greenyellow %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.greenyellow

#grey

modules = c("grey")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.grey <- read_delim("CytoscapeInputglom-edges0-grey.txt")
glom.edge.grey %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.grey
glom.edge.grey %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.grey

#grey60

modules = c("grey60")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.grey60 <- read_delim("CytoscapeInputglom-edges0-grey60.txt")
glom.edge.grey60 %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.grey60
glom.edge.grey60 %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.grey60

#lightcyan

modules = c("lightcyan")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.lightcyan <- read_delim("CytoscapeInputglom-edges0-lightcyan.txt")
glom.edge.lightcyan %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.lightcyan
glom.edge.lightcyan %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.lightcyan

#lightgreen

modules = c("lightgreen")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.lightgreen <- read_delim("CytoscapeInputglom-edges0-lightgreen.txt")
glom.edge.lightgreen %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.lightgreen
glom.edge.lightgreen %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.lightgreen

#lightyellow

modules = c("lightyellow")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.lightyellow <- read_delim("CytoscapeInputglom-edges0-lightyellow.txt")
glom.edge.lightyellow %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.lightyellow
glom.edge.lightyellow %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.lightyellow

#magenta

modules = c("magenta")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.magenta <- read_delim("CytoscapeInputglom-edges0-magenta.txt")
glom.edge.magenta %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.magenta
glom.edge.magenta %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.magenta

#midnightblue

modules = c("midnightblue")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.midnightblue <- read_delim("CytoscapeInputglom-edges0-midnightblue.txt")
glom.edge.midnightblue %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.midnightblue
glom.edge.midnightblue %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.midnightblue

#orange

modules = c("orange")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.orange <- read_delim("CytoscapeInputglom-edges0-orange.txt")
glom.edge.orange %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.orange
glom.edge.orange %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.orange

#paleturquoise

modules = c("paleturquoise")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.paleturquoise <- read_delim("CytoscapeInputglom-edges0-paleturquoise.txt")
glom.edge.paleturquoise %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.paleturquoise
glom.edge.paleturquoise %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.paleturquoise

#pink

modules = c("pink")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.pink <- read_delim("CytoscapeInputglom-edges0-pink.txt")
glom.edge.pink %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.pink
glom.edge.pink %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.pink

#purple

modules = c("purple")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);


glom.edge.purple <- read_delim("CytoscapeInputglom-edges0-purple.txt")
glom.edge.purple %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.purple
glom.edge.purple %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.purple

#red

modules = c("red")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.red <- read_delim("CytoscapeInputglom-edges0-red.txt")
glom.edge.red %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.red
glom.edge.red %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.red

#royalblue

modules = c("royalblue")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.royalblue <- read_delim("CytoscapeInputglom-edges0-royalblue.txt")
glom.edge.royalblue %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.royalblue
glom.edge.royalblue %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.royalblue

#saddlebrown

modules = c("saddlebrown")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.saddlebrown <- read_delim("CytoscapeInputglom-edges0-saddlebrown.txt")
glom.edge.saddlebrown %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.saddlebrown
glom.edge.saddlebrown %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.saddlebrown

#salmon

modules = c("salmon")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.salmon <- read_delim("CytoscapeInputglom-edges0-salmon.txt")
glom.edge.salmon %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.salmon
glom.edge.salmon %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.salmon

#skyblue

modules = c("skyblue")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.skyblue <- read_delim("CytoscapeInputglom-edges0-skyblue.txt")
glom.edge.skyblue %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.skyblue
glom.edge.skyblue %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.skyblue

#steelblue

modules = c("steelblue")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.steelblue <- read_delim("CytoscapeInputglom-edges0-steelblue.txt")
glom.edge.steelblue %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.steelblue
glom.edge.steelblue %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.steelblue

#tan

modules = c("tan")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.tan <- read_delim("CytoscapeInputglom-edges0-tan.txt")
glom.edge.tan %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.tan
glom.edge.tan %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.tan

#turquoise

modules = c("turquoise")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.turquoise <- read_delim("CytoscapeInputglom-edges0-turquoise.txt")
glom.edge.turquoise %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.turquoise
glom.edge.turquoise %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.turquoise

#violet

modules = c("violet")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.violet <- read_delim("CytoscapeInputglom-edges0-violet.txt")
glom.edge.violet %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.violet
glom.edge.violet %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.violet

#yellow

modules = c("yellow")
probes = colnames(wgcnaExpr) 
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputglom-edges0-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInputglom-nodes0-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

glom.edge.yellow <- read_delim("CytoscapeInputglom-edges0-yellow.txt")
glom.edge.yellow %>% dplyr::filter(fromNode %in% glom.degs$symbol) -> glom.edge.yellow
glom.edge.yellow %>% dplyr::filter(toNode %in% glom.degs$symbol) -> glom.edge.yellow


rbind(glom.edge.blue, glom.edge.brown, glom.edge.cyan, glom.edge.darkgreen,
      glom.edge.darkturquoise, glom.edge.green, glom.edge.greenyellow, glom.edge.lightcyan,
      glom.edge.lightyellow, glom.edge.magenta, glom.edge.midnightblue, glom.edge.orange,
      glom.edge.purple, glom.edge.red, glom.edge.steelblue, glom.edge.tan,
      glom.edge.turquoise, glom.edge.yellow) -> wgcna.glom.network


write_delim(wgcna.glom.network, "wgcna_glom_network.txt")
writexl::write_xlsx(wgcna.glom.network, "wgcna_glom_network.xlsx")


#Extract HUB-genes
clr = sort(unique(moduleColors))
hubs = character()

for (i in 1:length(clr)) {
  dt = geneModuleMembership[,i]
  kk = order(dt, decreasing = T)
  nm = rownames(geneModuleMembership)[kk]
  top10 = nm[1:10] #chose a number of top HUB-genes to extract from each module
  hubs = cbind(hubs,top10)
}

hubs = as.data.frame(hubs)
colnames(hubs) = clr

write.table(hubs, "Top_10_ModuleMembership_genes.txt", sep = "\t", quote = F)

hub = chooseTopHubInEachModule(wgcnaExpr, moduleColors, omitColors = "grey",
                               power = 1, type = "unsigned")

write.table(hub,"Module_top1_HUB-genes.txt", quote = F, col.names = F)

