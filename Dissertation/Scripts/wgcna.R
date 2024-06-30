rm(list = ls())
getwd()

library(WGCNA)
library(GEOquery)
library(tidyverse)

library(gridExtra)

allowWGCNAThreads()      

#Fetch Data for Type 2 diabetes gene matrix
data <- read.delim('diabetesrealcounts.txt',header= T )


#changing column name
colnames(data)
names(data) <- gsub('\\..*', '', names(data))

#checking for outliers
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK


#detecting heiharchial clustering
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

#detecting outliers in PCA
pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

colnames(data)

#modifying data column sequence

HCData <- data[, c(1:7 ,14:23, 34:42, 45:71, 74:78, 82:93,96)]
DMData <- data[,c(8:13,24:33,43,44,72,73,79:81,94,95)]
rt = cbind(HCData,DMData)
HCNum = ncol(HCData)
DMNum = ncol(DMData)





# Network Construction 
# choosing a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

data<-t(data)
rt<- t(rt)

# Call the network topology analysis function--------

# picking threshold value
sft <- pickSoftThreshold(rt, powerVector = power, verbose = 5)


sft.data <- sft$fitIndices

#visualizaton to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
 
#chose 6 as threshold value for my DM data

# memory estimate w.r.t blocksize
net = blockwiseModules(
  rt,
  power = 6,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  verbose = 3
)

#module eigengenes
MEs <- net$MEs


#get number of genes for each module
table(net$colors)
mergedColors = labels2colors(net$colors)
table(mergedColors)

# plotting color dendogram
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels =net$colors
moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]];
save(MEs,moduleLabels,moduleColors,geneTree,
     file = 'DMnetworkconstruction-auto.RData')

#Relating modules to traits


Type = c(rep('HC',HCNum),rep('DM',DMNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("HC","DM")
design = as.data.frame(design)
moduleColors <- labels2colors(net$colors)


#finding p value and correlation value of eignefactors
MEs0=moduleEigengenes(rt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
nSamples = nrow(rt)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)


#heatmap construction 
sizeGrWindow(10,8)
par(mar = c(2, 5.5, 2, 2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.43,
               cex.lab=0.5,
               cex.main=1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

modNames = substring(names(MEs), 3)

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.
geneModuleMembership = as.data.frame(cor(rt, MEs, use = "p"))
write.table(geneModuleMembership,"geneModuleMembershipDM1.csv",row.names=TRUE,col.names=TRUE,sep=",")

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))



#work in progress
names(geneModuleMembership) = paste("MM", modNames, sep="")

names(MMPvalue) = paste("p.MM", modNames, sep="")
DM=as.data.frame(design[,2])
names(DM) = "DM"
geneTraitSignificance = as.data.frame(cor(rt, DM, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));


names(geneTraitSignificance) = paste("GS.", names(DM), sep="");
write.table(geneTraitSignificance,'genetraitsignificancedm1.csv',row.names = TRUE, col.names = TRUE, sep = ',' )
names(GSPvalue) = paste("p.GS.", names(DM), sep="")
write.table(GSPvalue,"GSPvalue.csv",row.names=TRUE,col.names=TRUE,sep=",")
module = "green"
column = match(module, modNames);

moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(DMrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DM",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "green"
probes = colnames(rt) 
inModule = (moduleColors==module);
intModules = c("green")
for (module in intModules)
{
  modGenes = (moduleColors==module)
  modLLIDs = probes[modGenes];
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}


nGenes = ncol(rt)
nSamples = nrow(rt)
dissTOM <- 1-TOMsimilarityFromExpr(rt, power = 6);

nSelect = 500
set.seed(10)
select <- sample(nGenes, size = nSelect);
selectTOM <- dissTOM[select, select];
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select]


sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


MEs <- moduleEigengenes(rt, moduleColors)$eigengenes
DM=as.data.frame(design[,2]);
names(DM) = "DM"
MET <- orderMEs(cbind(MEs, DM))


sizeGrWindow(7,7);
par(cex = 0.9)
par(mar=c(4,5,1,1))
plotEigengeneNetworks(MET, "",
                      marDendro = c(0,4,1,2),
                      marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, xLabelsAngle
                      = 90)





