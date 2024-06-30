rm(list = ls())
options(stringsAsFactors = F)
getwd() 
 
library(ggplot2)
library(ggridges)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(tidyverse)
library(biomaRt)
data1=read.table("DIFF_all.DM.txt",header = T,sep="\t",check.names = FALSE)
names(data1)[1] <- 'Gene symbol'
#converting gene symbol to gene ID
listEnsembl()
ensembl <- useEnsembl(biomart = 'genes')
datasets <- listDatasets(ensembl)

ensembl.con <- useMart('ensembl',dataset = 'hsapiens_gene_ensembl')


attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

new_data <- getBM(attributes = c(	'external_gene_name','entrezgene_id'),
                  filters = 'external_gene_name',
                  values = data1$`Gene symbol`,
                  mart = ensembl.con)

#merging the entrez gene id to my data


data1 <- inner_join(data1, new_data, by = c("Gene symbol" = "external_gene_name"))
data1 <- data1[-c(1)]

#creating gene list
data1 <-na.omit(data1)
colnames(data1)
data1 <- data1[order(-data1$logFC),]


gene_list <- data1$logFC
names(gene_list) <- data1$entrezgene_id

gene_list


#gsea  GO analysis
gse <- gseGO(gene_list,ont = 'BP',
             keyType = 'ENTREZID',
             OrgDb = 'org.Hs.eg.db',
             eps = 1e-300)

gse1<- as.data.frame(gse)
write.table(gse1,file = 'gogsemf.txt',sep = '/t')
#gsea kegg analysis
KEGG_gseresult <- gseKEGG(gene_list, minGSSize = 10, maxGSSize = 1000)
kegg <- as.data.frame(KEGG_gseresult)

write.csv(kegg,'gseakeggt2d.csv')
#visualisation


par(mar=c(1,1,1,1))

png("gseaplot2_output.png", width=1200, height=1200)
ridgeplot(KEGG_gseresult,showCategory=10,label_format=2)
dev.off()

png("gseaplot2_output.png", width=1200, height=1200)
gseaplot2(KEGG_gseresult,1:10,pvalue_table = TRUE)
dev.off()



gseaplot(KEGG_gseresult,geneSetID = 2)
