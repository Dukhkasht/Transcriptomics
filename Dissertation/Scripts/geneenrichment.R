library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
gene_list <- c('SKIC3', 
               'CSRNP3', 
               'OPA1', 
               'C11orf58', 
               'EIF4G2', 
               'MATR3', 
               'TUSC3', 
               'PCLO', 
               'BTG3', 
               'SNAP25', 
               'CLCN4', 
               'SALL2', 
               'SLC17A6', 
               'NKX2-2', 
               'ISL1'
)
data1 <- data.frame(DM = gene_list)
GO_database <- 'org.Hs.eg.db'
gene2 <- bitr(data1$DM, fromType = 'SYMBOL', toType = 'GENENAME', OrgDb = GO_database)
GO <- enrichGO(gene2$GENENAME,
               OrgDb = GO_database,
               keyType = "GENENAME",
               ont = "ALL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2,
               readable = TRUE)
BP <- GO[GO$ONTOLOGY == 'BP', ]
CC <- GO[GO$ONTOLOGY == 'CC', ]
MF <- GO[GO$ONTOLOGY == 'MF', ]

write.table(as.data.frame(BP), 'go.DS.BP.txt', sep = '\t', row.names = T, quote = FALSE)
write.table(as.data.frame(CC), 'go.DS.CC.txt', sep = '\t', row.names = T, quote = FALSE)
write.table(as.data.frame(MF), 'go.DS.MF.txt', sep = '\t', row.names = T, quote = FALSE)


gene1 <- bitr(data1$DM,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
dotplot(GO, split="ONTOLOGY",showCategory=7,font.size=8,label_format=60)+facet_grid(ONTOLOGY~., scale="free")
              
KEGG<-enrichKEGG(gene1$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
keggdata <- as.data.frame(KEGG)
barplot(KEGG,showCategory=40,font.size=8,label_format=60, title = 'kegg enrichment for Common genes in T2D and Male infertility')
save(KEGG,file = 'KEGG DS.Rdata')
write.csv(keggdata,'keggenrichment.csv')
