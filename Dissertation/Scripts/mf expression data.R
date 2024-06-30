rm(list =ls())
getwd()
library(GEOquery)
library(affy)
#getting series matrix file from GEO
eSet <- getGEO("GSE6872",
                destdir = '.',
                getGPL = F) 

#creating expression file
gse6872 <- as.data.frame(exprs(eSet[[1]]))

#getting file using affy library
untar('GSE6872_RAW.tar',exdir = 'data/')

raw.data <- ReadAffy(celfile.path = "data/")
eset <- rma(raw.data)
gse6872 <- as.data.frame(exprs(eset))





#getting phenotype data
gse <- getGEO('gse6872',GSEMatrix  = TRUE)
pheno <- pData(phenoData(gse[[1]]))


#annotating using soft file
GPL=getGEO(filename = 'gse6872_family.soft.gz') 
gpl=GPL@gpls[[1]]@dataTable@table
colnames(gpl)

ids <- gpl[,c(1,11)]
rownames(ids) <- ids$ID
colnames(ids) = c("probe_id" ,"symbol")
gse6872$probe_id=rownames(gse6872)
gse6872=merge(gse6872,ids,by.x="probe_id", by.y="probe_id") 
gse6872=gse6872[!duplicated(gse6872$symbol),]  
rownames(gse6872)=gse6872$symbol
gse6872=gse6872[,-c(1)]
colnames(gse6872)
save(gse6872,file = 'Expression6872.Rdata')


### using annotation files
data1<- read.csv('Id35.csv',sep=',',header=T)
gse6872$ID = rownames(gse6872)

gse6872= merge(gse6872,data1,by.x='ID',by.y='ID')
colnames(gse6872)
gse6872 = gse6872[!duplicated(gse6872$Gene.symbol),]
rownames(gse6872) = gse6872$Gene.symbol
gse6872 = gse6872[,-c(1)]
save(gse6872,file = 'Expression45887.Rdata')

####

