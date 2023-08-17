
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
####logFC为2
m=ggplot(data=allDiff, aes(x=logFC, y =-log10(FDR))) +
  #将数据分成四象限，对每个象限的数据进行颜色和大小的定义
  geom_point(data=subset(data,abs(data$logFC) <= 2),aes(size=abs(logFC)),color="gray",alpha=0.5) +
  geom_point(data=subset(data,data$FDR >= 0.05&abs(data$logFC) > 2),aes(size=abs(logFC)),color="gray",alpha=0.5) +
  geom_point(data=subset(data,data$FDR<0.05 & data$logFC > 2),aes(size=abs(logFC)),color="red",alpha=0.5) +
  geom_point(data=subset(data,data$FDR<0.05 & data$logFC < -2),aes(size=abs(logFC)),color="darkgreen",alpha=0.5) +
  ## 添加分割线（横向+竖向）
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(2,-2),lty=4,lwd=0.6,alpha=0.8)+
  #设置背景为空白
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  ## 修改坐标轴
  labs(x="Log2FC",y="-Log10 (FDR)")+
  ## 去掉图注
  theme(legend.position='none')
## 添加文本信息,此处为添加logFC>5的基因名
#geom_text_repel(data=subset(data, abs(logFC) > 5), aes(label=id),color="black",alpha = 0.8)
print(m) 
ggsave(m,file="volcano.pdf",width = 8,height = 6)#保存为pdf格式

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
         ) #不显示列名
dev.off()

#add entrezIDs

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

library("org.Hs.eg.db")          #引用包
rtdata=read.table("diffSig.txt",sep="\t",check.names=F,header=T)    #读取文件
#rownames to column
rtdata=rtdata[,c(1,2)]
write.table(rtdata,file="geneID.txt",sep="\t",quote=F,row.names=F)    #输出结果

genes=as.vector(rtdata[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #找出基因对应的id
entrezIDs <- as.character(entrezIDs)
out=cbind(rtdata,entrezID=entrezIDs)
write.table(out,file="geneEntrezID.txt",sep="\t",quote=F,row.names=F)    #输出结果

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
rt=rt[is.na(rt[,"entrezID"])==F,]                                 #去除基因id为NA的基因
gene=rt$entrezID

#GO富集分析
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果

#柱状图
pdf(file="barplot-go.pdf",width = 10,height = 10)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#气泡图
pdf(file="bubble-go.pdf",width = 10,height = 10)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#kegg富集分析
rt=rt[is.na(rt[,"entrezID"])==F,]                             #去除基因id为NA的基因
gene=rt$entrezID

#富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #富集分析
write.table(kk,file="KEGGId.txt",sep="\t",quote=F,row.names = F)                          #保存富集结果

#柱状图
pdf(file="barplot-kegg.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()

#气泡图
pdf(file="bubble-kegg.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 30)
dev.off()

#GSEA
rt=read.table("geneEntrezID.txt",sep="\t",header=T,check.names=F)   
rt=rt[is.na(rt[,"entrezID"])==F,]  
head(rt)
# 在对数据按照logFC从大到小排序，这里还可以使用dplyr的arrange

gsea_data_order <- rt[order(as.vector(rt$logFC), decreasing = T),]
gsea_data_order <- gsea_data_order[is.na(gsea_data_order[,"entrezID"])==F,]
head(gsea_data_order)
#整理输入数据
gene_logFC <- gsea_data_order$logFC
gene_id <- gsea_data_order$entrezID
names(gene_logFC) <- gene_id
head(gene_logFC)

#GSEA-GO富集分析
#ont可改为BP、CC、MF或ALL
go <- gseGO(gene_logFC,OrgDb= org.Hs.eg.db,ont= "ALL",minGSSize= 100,maxGSSize= 500,pvalueCutoff = 0.05,verbose= FALSE)
write.table(go,file="GSEA-ALL.txt",sep="\t",quote=F,row.names = F)
##对GSEA结果进行可视化：
go_order<-go[order(go$NES,decreasing=T)] #按enrichmentScore或NES降序排列
gseaplot2(go, row.names(go_order)[1:5],base_size = 10) #前5个GO
gseaplot2(go, row.names(go_order)[1:5],base_size = 10,pvalue_table = TRUE) #前5个GO,加p值
gseaplot2(go, row.names(go_order)[3],base_size = 10,title = go$Description[3]) #第3个GO
gseaplot2(go, "GO:0000278",base_size = 10,color='steelblue') #换颜色
#将结果画成山脊图
ridgeplot(go,showCategory = 30, fill = "pvalue")

#GSEA-KEGG富集分析
kk<-gseKEGG(gene_logFC, organism = "hsa",pvalueCutoff = 0.05,
            nPerm = 1000, minGSSize = 10, maxGSSize = 500,
            verbose = TRUE, seed = FALSE, by = "fgsea")
kk_order<-kk[order(kk$enrichmentScore,decreasing=T)]
# 写出富集分析结果
write.table(kk_order,file="GSEA-KEGG.txt",sep="\t",quote=F,row.names = F)

##对GSEA结果进行可视化：
gseaplot2(kk, row.names(kk_order)[1:5],base_size = 10)
gseaplot2(kk, row.names(kk_order)[3],base_size = 10,title = kk$Description[3]) #第3个GO,改俩地方
#将结果画成山脊图
ridgeplot(kk,fill = "pvalue")



diffFile="geneID.txt"            #差异分析的结果文件
gmtFile="immunesigdb.gmt"     #免疫基因集文件

#读取文件,并对输入文件进行整理
rt=read.table(diffFile, header=T, sep="\t", check.names=F)
rt=rt[order(rt$logFC, decreasing=T),]
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])

#读入基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff=1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

#绘制实验组富集的图形
termNum=5     #展示前5个基因集
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]       #需要展示的基因集名称
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in treat group")
  pdf(file="GSEA.treat.pdf", width=8, height=6)
  print(gseaplot)
  dev.off()
}

#绘制对照组富集的图形
termNum=5     #展示前5个基因集
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]       #需要展示的基因集名称
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in control group")
  pdf(file="GSEA.con.pdf", width=8, height=6)
  print(gseaplot)
  dev.off()
}
 
