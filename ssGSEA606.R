if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
#BiocManager::install("rtracklayer")

rm(list=ls())
setwd("/Users/sunx/Documents/ovca-immune/ucsc/NMF606/ssGSEA/")
library(dplyr)
library(tidyverse)

gz=gzfile('/Users/sunx/Documents/ovca-immune/ucsc/NMF606/ssGSEA/PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level.txt.gz','rt')
dataOrg=read.table(gz,sep = "\t",header = T,check.names = F)

# import readr
#TCGASubtypeOrg <- read_delim("~/Documents/ovca-immune/ucsc/explore/subtype/TCGASubtype.20170308.tsv.gz", 
#                            delim = "\t", escape_double = FALSE, 
#                            trim_ws = TRUE)
data=dataOrg
data[1:2,1:2]
colnames(data)[1]='sample'
names=colnames(data)
types=read.table('TCGA_GTEX_category.txt',sep = '\t',header = T,check.names = F)
#ovary=types %>% filter(str_detect(TCGA_GTEX_main_category,'Ovary'))
ovarian=types %>% filter(str_detect(TCGA_GTEX_main_category,'Ovarian'))
#ovaryID=c(ovary[,1])
#ovarydata=select(data,ovaryID)

ovarianID=ovarian[,1]
ovarianID=c('sample',ovarianID)
interID=intersect(names,ovarianID)
ssGSEAOvca=select(data,interID)
colnames(ssGSEAOvca)[1]='id'
riskdata=read.table('risk23.txt',sep = '\t',header = T,check.names = F)
inter=intersect(interID,riskdata$id)
rownames(ssGSEAOvca)=ssGSEAOvca[,1]
ssGSEAOvca=ssGSEAOvca[,-1]
ssGSEAOvca=as.data.frame(t(ssGSEAOvca))
ssGSEAOvca=cbind(id=rownames(ssGSEAOvca),ssGSEAOvca)
ssGSEARiskOvca <- inner_join(ssGSEAOvca, riskdata,
                      by = "id") 

write.table(ssGSEARiskOvca, file="ssGSEARiskOvca.txt", sep="\t", quote=F, row.names=F, col.names=T)

#重新读入数据
riskssGSEA=read.table('ssGSEARiskOvca.txt',header = T, sep = '\t',check.names = F)
tail(colnames(riskssGSEA))
riskssGSEA1=riskssGSEA[,c(1:1387,ncol(riskssGSEA))]
riskssGSEAb=riskssGSEA[,c(1:1387,(ncol(riskssGSEA)-1),ncol(riskssGSEA))]
tail(colnames(riskssGSEA1))
tail(colnames(riskssGSEAb))
riskssGSEA1=riskssGSEA1[order(riskssGSEA1$risk),]
table(riskssGSEA1$risk)
group=riskssGSEA1[,c(1,ncol(riskssGSEA1))]
rownames(group)=group[,1]
group=group$risk
#group=as.factor(group)

riskssGSEA2=riskssGSEA1
rownames(riskssGSEA2)=riskssGSEA2$id
riskssGSEA2=riskssGSEA2[,-c(1,ncol(riskssGSEA2))]
tail(colnames(riskssGSEA2))
riskssGSEA2=t(t(riskssGSEA2))

##################################################
#limma找差异
library(limma)

riskssGSEAa=riskssGSEA1
rownames(riskssGSEAa)=riskssGSEAa$id
riskssGSEAa=riskssGSEAa[,-c(1,ncol(riskssGSEAa))]
tail(colnames(riskssGSEAa))
riskssGSEAa=t(riskssGSEAa)

# 设定分组
condition <- factor(group)
table(condition)
design1=model.matrix(~factor(condition))
fit1=lmFit(riskssGSEAa,design1)
fit1=eBayes(fit1)
options(digits = 4)
b2 <- topTable(fit1,coef=2,adjust='BH',n=Inf)
#b3 <- topTable(fit1,coef=3,adjust='BH')
#b4 <- topTable(fit1,coef=4,adjust='BH',n=Inf)
length(which(b2$adj.P.Val < 0.05))
length(which(b2$adj.P.Val < 0.02))
top72=b2[b2$adj.P.Val < 0.05,]
top17=b2[b2$adj.P.Val < 0.02,]
top17list=rownames(top17)
write.table(top17list,file="top17list.txt",sep="\t",quote=F,col.names=F,row.names = F) 


#分组聚类！！！
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

rownames(riskssGSEAb)=riskssGSEAb$id
riskssGSEAb=riskssGSEAb[,-c(1,ncol(riskssGSEAb))]
tail(colnames(riskssGSEAb))
riskssGSEAb=riskssGSEAb[order(riskssGSEAb$riskScore,decreasing = T),] #排序很重要
riskssGSEAb=t(riskssGSEAb)

top17list1=c('riskScore',top17list)
ssGSEAb.mat=riskssGSEAb[rownames(riskssGSEAb) %in% top17list1,]

Group=group
table(Group)
col_fun = colorRamp2(c(-0.3, 0, 0.3),c("#2fa1dd", "white", "red"))
#col = colorRamp2(breaks = seq(min(rt1), max(rt1), length = 3), colors = c("#2fa1dd", "white", "#f87669"))

top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                       labels = c("low","high"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

pdf(file="ssGSEA17.03.pdf",width = 15,height = 6)
Heatmap(ssGSEAb.mat,name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = Group,
        show_heatmap_legend = T,
        border = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL)
dev.off()



###################################################
#筛选有差异的，用LPS做t检验
library(LPS)
coeff <- LPS.coeff(data=riskssGSEA2, response=group) #from t test value
coeff1=as.data.frame(coeff)
sigcoeff005=coeff[as.numeric(as.vector(coeff1$p.value))<0.05,] 
sigcoeff001=coeff[as.numeric(as.vector(coeff1$p.value))<0.01,] 
sigcoeff005a=cbind(id=rownames(sigcoeff005),sigcoeff005)
sigcoeff001a=cbind(id=rownames(sigcoeff001),sigcoeff001)
write.table(sigcoeff005a,file="ssGSEA005.Sig.txt",sep="\t",row.names=F,quote=F)
write.table(sigcoeff001a,file="ssGSEA001.Sig.txt",sep="\t",row.names=F,quote=F)

m <- LPS(data=riskssGSEA2, coeff=coeff, response=group, k=30)

pdf('ssGSEA001.pdf',width = 10,height = 12)
predict(m, riskssGSEA2, type="probability",threshold = 0.8, plot=TRUE)
dev.off()

riskssGSEA3=riskssGSEA2[,colnames(riskssGSEA2) %in% rownames(sigcoeff005)]
riskssGSEA3=t(riskssGSEA3)

riskssGSEA4=riskssGSEA2[,colnames(riskssGSEA2) %in% rownames(sigcoeff001)]
riskssGSEA4=t(riskssGSEA4)
library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(-0.3, 0, 0.3),c("#2fa1dd", "white", "red"))
#col = colorRamp2(breaks = seq(min(rt1), max(rt1), length = 3), colors = c("#2fa1dd", "white", "#f87669"))

top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c( "#f87669","#2fa1dd")),
                       labels = c("high","low"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

pdf(file="ssGSEArisk005.pdf",width = 15,height = 6)
Heatmap(riskssGSEA3,name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = group,
        show_heatmap_legend = T,
        border = F,
        show_column_names = F,
        show_row_names = F,
        column_title = NULL)
dev.off()

pdf(file="ssGSEArisk001.pdf",width = 15,height = 6)
Heatmap(riskssGSEA4,name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = group,
        show_heatmap_legend = T,
        border = F,
        show_column_names = F,
        show_row_names = F,
        column_title = NULL)
dev.off()

library(vioplot) 

#对样品分组，用0.01的121个来做
conNum=122
treatNum=127
rt=riskssGSEA4
conData=rt[,c(1:122)]
treatData=rt[,c(123:ncol(rt))]
rt=t(rt)
#绘制小提琴图
outTab=data.frame()
pdf(file="vioplot.pdf", width=80, height=60)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,(3*ncol(rt)-3)),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=10,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，对照组用蓝色表示，实验组用红色表示
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Con", "Treat"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))
text(seq(1,(3*ncol(rt)-2),3),-0.075,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

#输出0.01免疫细胞和p值表格文件
write.table(outTab,file="immuneDiff.txt",sep="\t",row.names=F,quote=F)


#有差异的重新做图
#对样品分组
immuDiff=read.table('immuneDiff.txt',header = T, sep = '\t',check.names = F)
conNum=122
treatNum=127
riskssGSEA5=riskssGSEA4
rtdif=riskssGSEA5[rownames(riskssGSEA5) %in% immuDiff$Cell,]
rt=rtdif
conData=rt[,c(1:122)]
treatData=rt[,c(123:ncol(rt))]
rt=t(rt)
#绘制小提琴图
outTab=data.frame()
pdf(file="vioplot-immudif.pdf", width=15, height=10)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,(3*ncol(rt)-3)),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=8,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，对照组用蓝色表示，实验组用红色表示
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'red')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'blue')
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("high", "low"),
       lwd=3,bty="n",cex=1,
       col=c("red","blue"))
text(seq(1,(3*ncol(rt)-2),3),-0.075,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

# 分组信息
highriskhigh=outTab[c(2,3,4,6,13,19),]
lowriskhigh=outTab[-c(2,3,4,6,13,19),]
write.table(highriskhigh,file="highriskhigh.txt",sep="\t",row.names=F,quote=F)
write.table(lowriskhigh,file="lowriskhigh.txt",sep="\t",row.names=F,quote=F)

#ssGSEA 23 cell heatmap
col_fun = colorRamp2(c(-0.3, 0, 0.3),c("#2fa1dd", "white", "red"))
#col = colorRamp2(breaks = seq(min(rt1), max(rt1), length = 3), colors = c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c( "#f87669","#2fa1dd")),
                       labels = c("high","low"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

pdf(file="ssGSEA23risk.pdf",width = 15,height = 6)
Heatmap(rtdif,name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = group,
        show_heatmap_legend = T,
        border = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL)
dev.off()


#用0.05的riskssGSEA3
conNum=122
treatNum=127
rt=riskssGSEA3
conData=rt[,c(1:122)]
treatData=rt[,c(123:ncol(rt))]
rt=t(rt)
#绘制小提琴图
outTab=data.frame()
pdf(file="vioplot.pdf", width=80, height=60)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,(3*ncol(rt)-3)),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=10,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，对照组用蓝色表示，实验组用红色表示
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Con", "Treat"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))
text(seq(1,(3*ncol(rt)-2),3),-0.075,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

#输出0.05免疫细胞和p值表格文件
write.table(outTab,file="immuneDiff005.txt",sep="\t",row.names=F,quote=F)


#对0.05有差异的重新做图
#对样品分组
immuDiff005=read.table('immuneDiff005.txt',header = T, sep = '\t',check.names = F)
conNum=122
treatNum=127
riskssGSEA6=riskssGSEA3
rtdif1=riskssGSEA6[rownames(riskssGSEA6) %in% immuDiff005$Cell,]
rt=rtdif1
conData=rt[,c(1:122)]
treatData=rt[,c(123:ncol(rt))]
rt=t(rt)
#绘制小提琴图
outTab=data.frame()
pdf(file="vioplot-immudif005.pdf", width=15, height=10)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,(3*ncol(rt)-3)),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=8,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，对照组用蓝色表示，实验组用红色表示
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'red')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'blue')
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("high", "low"),
       lwd=3,bty="n",cex=1,
       col=c("red","blue"))
text(seq(1,(3*ncol(rt)-2),3),-0.075,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()






#WGCNA
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

expro=read.csv('immurisk.txt',sep = '\t',row.names = 1)
colnames(expro)
immue68=expro[,1:68]
riskfactor=expro[,c(69,71:94)]

#相关性分析
outTab=data.frame()
for(cell in colnames(immue68)){
  for(gene in colnames(riskfactor)){
    x=as.numeric(immue68[,cell])
    y=as.numeric(riskfactor[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
  }
}

#绘制相关性热图
outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=14, height=12)
ggplot(outTab, aes(Gene, Immune)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    #去掉背景
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   #x轴字体
        axis.text.y = element_text(size = 8, face = "bold")) +       #y轴字体
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #设置图例
  scale_x_discrete(position = "bottom")      #X轴名称显示位置
dev.off()





