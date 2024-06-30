rm(list = ls())
options(stringsAsFactors = F)

library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
library(grid)
rt=read.table("diabetesrealcounts.txt",sep="\t",header=T,check.names=F,row.names = 1)

data=rt


#for male infertitlity data
HCData <- data[, c(1:13,27, 31:34, 41)]

MFData <- data[, c(14:26,28:30, 35:40)]

#for diabetes data
HCData <- data[, c(1:7,14:23, 34:42, 45:71,74:78,82:93,96)]

DMData <- data[, c(8:13,24:33, 43:44,72:73,79:81,94,95)]


rt=cbind(HCData,DMData)
HCNum=ncol(HCData)
DMNum=ncol(DMData)
Type=c(rep("HC",HCNum),rep("DM",DMNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("HC","DM")
fit <- lDMit(rt,design)
cont.matrix<-makeContrasts(DM-HC,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

DiffDMDM=topTable(fit2,adjust='fdr',number=length(rownames(data)))
write.table(DiffDMDM,file="DiffDM_all.DM.txt",sep="\t")



DiffDM <-  read.table('Diff_all.DM.txt',sep = '\t',row.names = 1,header = T)
#---------------------------------#
#visualisation using heatmap and volcano plot

DiffDM=DiffDM[order(as.numeric(as.vector(DiffDM$logFC))),]
DiffDMGene=as.vector(rownames(DiffDM))
DiffDMLength=length(DiffDMGene)
afGene=c()
if(DiffDMLength>(100)){
  afGene=DiffDMGene[c(1:50,(DiffDMLength-50+1):DiffDMLength)]
}else{
  afGene=DiffDMGene
}
afExp=rt[afGene,]
n=t(scale(t(afExp)))
n[n>2]=2
n[n< -2]= -2
Type=c(rep("HC",HCNum),rep("DM",DMNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
anncolor=list(Type=c(DM=pal_npg()(1),HC=pal_npg()(2)[2]))
##pdf(file="DiffDMDM_heatmap.pdf",height=7,width=8)
pheatmap(n,                                                                      
         annotation=Type,                                                            
         color = colorRampPalette(c(pal_npg()(2)[2],"white", pal_npg()(1)))(50),     
         cluster_cols =F,                                                           
         show_colnames = F,                                                         
         scale="row", 
         fontsize = 8,
         fontsize_row=3,
         fontsize_col=8,
         annotation_colors=anncolor
)
#dev.off()
P=0.05
aflogFC=0.175
Significant=ifelse((DiffDM$P.Value<P & abs(DiffDM$logFC)>aflogFC), ifelse(DiffDM$logFC>aflogFC,"Up","Down"), "Not")
p = ggplot(DiffDM, aes(logFC, -log10(P.Value)))+
  geom_point(aes(col=Significant),size=3)+
  scale_color_manual(values=c(pal_npg()(2)[2], "#838B8B", pal_npg()(1)))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  geom_hline(aes(yintercept=-log10(P)), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=aflogFC), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=-aflogFC), colour="gray", linetype="twodash",size=1)
DiffDM$symbol=rownames(DiffDM)
png("DiffDM_vol.png",width=1200,height=1200)
p=p+theme_bw()
p+geom_point(size = 2, shape = 1,) +
  ggrepel::geom_label_repel(
    aes(label = DiffDM$symbol),
    color="black",
    label.size =0.1
  )
dev.off()



# --------------------- #
p <- pheatmap(n,
              annotation = Type,
              color = colorRampPalette(c("skyblue", "white", "red"))(50),
              cluster_cols = FALSE,
              show_colnames = FALSE,
              scale = "row", 
              fontsize = 8,
              fontsize_row = 5,
              fontsize_col = 8,
              annotation_colors = anncolor)
g <- p$gtable



  # Draw the original heatmap
png(filename = 'heatmapDM.png',width = 1000,height = 1000)
grid.draw(p$gtable)
dev.off()

