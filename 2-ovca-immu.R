
rm(list=ls())
setwd("/Users/sunx/Documents/ovca-immune/2pas-resis/")


library(limma)
library(edgeR)
# -------------------------------------------------------->>>>>>>>>>
# make
# -------------------------------------------------------->>>>>>>>>>

count_raw <- read.table(file = "/Users/sunx/Documents/ovca-immune/GSE163854_counts.txt",header = T,sep = "\t")
count_exp <- count_raw[,c(1,6,7,11,12,13)]
count_exp=as.matrix(count_exp)
rownames(count_exp)=count_exp[,1]
count_exp=count_exp[,2:ncol(count_exp)]
dimnames=list(rownames(count_exp),colnames(count_exp))
data=matrix(as.numeric(as.matrix(count_exp)),nrow=nrow(count_exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
colnames(data)
data=subset(data,select = c(3,4,5,1,2))
colnames(data)
group_info = c(rep("ctrl",3), rep("treated",2))
design <- model.matrix(~group_info)

y <- DGEList(counts=data,group=group_info)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("ctrl","treated"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

diffa=cbind(id=rownames(diff),diff)

write.table(diffa,file="all.txt",sep="\t",quote=F,col.names=T,row.names = F)
logFCcutoff=1
padj=0.05
diffSig = diff[(diff$FDR < padj & (diff$logFC>logFCcutoff | diff$logFC<(-logFCcutoff))),]
diffSiga=cbind(id=rownames(diffSig),diffSig)

write.table(diffSiga, file="diffSig.txt",sep="\t",quote=F,row.names = F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>logFCcutoff)),]
diffUpa=cbind(id=rownames(diffUp),diffUp)

write.table(diffUpa, file="up.txt",sep="\t",quote=F,row.names = F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-logFCcutoff))),]
diffDowna=cbind(id=rownames(diffDown),diffDown)
write.table(diffDowna, file="down.txt",sep="\t",quote=F,row.names = F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F) 
heatmapData <- newData[rownames(diffSig),]

#volcano
library(ggplot2)
library(ggrepel)
data = read.table("all.txt", header=TRUE)
####logFCΪ2
m=ggplot(data=allDiff, aes(x=logFC, y =-log10(FDR))) +
  #�����ݷֳ������ޣ���ÿ�����޵����ݽ�����ɫ�ʹ�С�Ķ���
  geom_point(data=subset(data,abs(data$logFC) <= 2),aes(size=abs(logFC)),color="gray",alpha=0.5) +
  geom_point(data=subset(data,data$FDR >= 0.05&abs(data$logFC) > 2),aes(size=abs(logFC)),color="gray",alpha=0.5) +
  geom_point(data=subset(data,data$FDR<0.05 & data$logFC > 2),aes(size=abs(logFC)),color="red",alpha=0.5) +
  geom_point(data=subset(data,data$FDR<0.05 & data$logFC < -2),aes(size=abs(logFC)),color="darkgreen",alpha=0.5) +
  ## ���ӷָ��ߣ�����+����
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(2,-2),lty=4,lwd=0.6,alpha=0.8)+
  #���ñ���Ϊ�հ�
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  ## �޸�������
  labs(x="Log2FC",y="-Log10 (FDR)")+
  ## ȥ��ͼע
  theme(legend.position='none')
## �����ı���Ϣ,�˴�Ϊ����logFC>5�Ļ�����
#geom_text_repel(data=subset(data, abs(logFC) > 5), aes(label=id),color="black",alpha = 0.8)
print(m) 
ggsave(m,file="volcano.pdf",width = 8,height = 6)#����Ϊpdf��ʽ

#heatmap
library(pheatmap)
hmExp=log10(heatmapData+0.001)
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf")
pheatmap(hmMat, 
         #color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =T,
         fontsize_row=3,
         fontsize_col=5,
         #show_rownames=F,
         #show_colnames=F
         ) #����ʾ����
dev.off()

#add entrezIDs

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

library("org.Hs.eg.db")          #���ð�
rtdata=read.table("diffSig.txt",sep="\t",check.names=F,header=T)    #��ȡ�ļ�
#rownames to column
rtdata=rtdata[,c(1,2)]
write.table(rtdata,file="geneID.txt",sep="\t",quote=F,row.names=F)    #������

genes=as.vector(rtdata[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #�ҳ������Ӧ��id
entrezIDs <- as.character(entrezIDs)
out=cbind(rtdata,entrezID=entrezIDs)
write.table(out,file="geneEntrezID.txt",sep="\t",quote=F,row.names=F)    #������

#go analysis
install.packages("colorspace")
install.packages("stringi")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("enrichplot")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

rt=read.table("geneEntrezID.txt",sep="\t",header=T,check.names=F)                   
rt=rt[is.na(rt[,"entrezID"])==F,]                                 #ȥ������idΪNA�Ļ���
gene=rt$entrezID

#GO��������
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 #���渻�����

#��״ͼ
pdf(file="barplot-go.pdf",width = 10,height = 10)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#����ͼ
pdf(file="bubble-go.pdf",width = 10,height = 10)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#kegg��������
rt=rt[is.na(rt[,"entrezID"])==F,]                             #ȥ������idΪNA�Ļ���
gene=rt$entrezID

#��������
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #��������
write.table(kk,file="KEGGId.txt",sep="\t",quote=F,row.names = F)                          #���渻�����

#��״ͼ
pdf(file="barplot-kegg.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()

#����ͼ
pdf(file="bubble-kegg.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 30)
dev.off()

#GSEA
rt=read.table("geneEntrezID.txt",sep="\t",header=T,check.names=F)   
rt=rt[is.na(rt[,"entrezID"])==F,]  
head(rt)
# �ڶ����ݰ���logFC�Ӵ�С�������ﻹ����ʹ��dplyr��arrange

gsea_data_order <- rt[order(as.vector(rt$logFC), decreasing = T),]
gsea_data_order <- gsea_data_order[is.na(gsea_data_order[,"entrezID"])==F,]
head(gsea_data_order)
#������������
gene_logFC <- gsea_data_order$logFC
gene_id <- gsea_data_order$entrezID
names(gene_logFC) <- gene_id
head(gene_logFC)

#GSEA-GO��������
#ont�ɸ�ΪBP��CC��MF��ALL
go <- gseGO(gene_logFC,OrgDb= org.Hs.eg.db,ont= "ALL",minGSSize= 100,maxGSSize= 500,pvalueCutoff = 0.05,verbose= FALSE)
write.table(go,file="GSEA-ALL.txt",sep="\t",quote=F,row.names = F)
##��GSEA������п��ӻ���
go_order<-go[order(go$NES,decreasing=T)] #��enrichmentScore��NES��������
gseaplot2(go, row.names(go_order)[1:5],base_size = 10) #ǰ5��GO
gseaplot2(go, row.names(go_order)[1:5],base_size = 10,pvalue_table = TRUE) #ǰ5��GO,��pֵ
gseaplot2(go, row.names(go_order)[3],base_size = 10,title = go$Description[3]) #��3��GO
gseaplot2(go, "GO:0000278",base_size = 10,color='steelblue') #����ɫ
#���������ɽ��ͼ
ridgeplot(go,showCategory = 30, fill = "pvalue")

#GSEA-KEGG��������
kk<-gseKEGG(gene_logFC, organism = "hsa",pvalueCutoff = 0.05,
            nPerm = 1000, minGSSize = 10, maxGSSize = 500,
            verbose = TRUE, seed = FALSE, by = "fgsea")
kk_order<-kk[order(kk$enrichmentScore,decreasing=T)]
# д�������������
write.table(kk_order,file="GSEA-KEGG.txt",sep="\t",quote=F,row.names = F)

##��GSEA������п��ӻ���
gseaplot2(kk, row.names(kk_order)[1:5],base_size = 10)
gseaplot2(kk, row.names(kk_order)[3],base_size = 10,title = kk$Description[3]) #��3��GO,�����ط�
#���������ɽ��ͼ
ridgeplot(kk,fill = "pvalue")



diffFile="geneID.txt"            #��������Ľ���ļ�
gmtFile="immunesigdb.gmt"     #���߻����ļ�

#��ȡ�ļ�,���������ļ���������
rt=read.table(diffFile, header=T, sep="\t", check.names=F)
rt=rt[order(rt$logFC, decreasing=T),]
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])

#��������ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff=1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

#����ʵ���鸻����ͼ��
termNum=5     #չʾǰ5������
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]       #��Ҫչʾ�Ļ�������
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in treat group")
  pdf(file="GSEA.treat.pdf", width=8, height=6)
  print(gseaplot)
  dev.off()
}

#���ƶ����鸻����ͼ��
termNum=5     #չʾǰ5������
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]       #��Ҫչʾ�Ļ�������
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in control group")
  pdf(file="GSEA.con.pdf", width=8, height=6)
  print(gseaplot)
  dev.off()
}
 