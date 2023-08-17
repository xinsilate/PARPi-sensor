
rm(list=ls())

setwd("/Users/sunx/Documents/ovca-immune/ucsc/NMF606/")


library(NMF)
#rt2=read.table("blue462.txt",sep="\t",header=T,check.names=F,row.names = 1)

rt1=read.table("gene606OS1_log2.txt",sep="\t",header=T,check.names=F,row.names = 1)
rt=rt1[,-c(1,2)]
rt=t(rt)
ranks <- 2:5
#运行很慢，最好用服务???
estim.nmf <- nmf(rt,ranks, nrun=50)
#选择合的rank
plot(estim.nmf)

#用rank2再次NMF
seed = 123
estim.nmf2 <- nmf(rt, rank = 2, nrun=50, seed = seed, method = "brunet")

#绘制分群??????
#设置颜色
jco <- c("pink","green",'blue')
index <- extractFeatures(estim.nmf2,"max") 
sig.order <- unlist(index)
NMF.Exp.rank <- rt[sig.order,]
NMF.Exp.rank <- na.omit(NMF.Exp.rank) #sig.order有时候会有缺失???
group <- predict(estim.nmf2) # 提出亚型
table(group)
write.table(group,file="group.txt",sep="\t",row.names=F,quote=F)

pdf("consensusmap.pdf",onefile = FALSE,
    width = 8,
    height =8)
consensusmap(estim.nmf2,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(NMF.Exp.rank)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2])))
dev.off()

pdf("basismap.pdf",onefile = FALSE,
    width = 8,
    height =8)
basismap(estim.nmf2,
         cexCol = 1.5,
         cexRow = 1,
         annColors=list(c("#2874C5","#EABF00","#C6524A")))
dev.off()


library(tidyverse)
groupid=read.table("group.txt",sep="\t",header=T,check.names=F)
expOSid=cbind(rt1,groupid=groupid$x)

library(survival)
library(survminer)
sfit <- survfit(Surv(OS.time, OS) ~ groupid,
                data = expOSid)

pdf("surv2.pdf",onefile = FALSE,
    width = 8,
    height =8)
ggsurvplot(sfit,pval = T,palette = "jco")
dev.off()

library(Rtsne)
tsne_out = Rtsne(t(rt),perplexity = 30)
pdat = data.frame(tsne_out$Y,factor(group))
colnames(pdat) = c("Y1","Y2","group")
head(pdat)

library(ggplot2)
library(paletteer)
pdf("tsne.pdf",onefile = FALSE,
    width = 8,
    height =8)
ggplot(pdat,aes(Y1,Y2))+
  geom_point(aes(Y1,Y2,fill = group),shape = 21,color = "black")+
  stat_ellipse(aes(color = group,fill = group),
               geom = "polygon",
               alpha = 0.3,
               linetype = 2)+
  scale_color_paletteer_d("RColorBrewer::Set3")+
  scale_fill_paletteer_d("RColorBrewer::Set3")+
  theme_classic()+
  theme(legend.position = "top")
dev.off()

#用rank4再次NMF
seed = 111
estim.nmf2 <- nmf(rt, rank = 4, nrun=50, seed = seed, method = "brunet")

#绘制分群
#设置颜色
jco <- c("deepskyblue3","lightgoldenrod2",'darkred','green4')
index <- extractFeatures(estim.nmf2,"max") 
sig.order <- unlist(index)
NMF.Exp.rank <- rt[sig.order,]
NMF.Exp.rank <- na.omit(NMF.Exp.rank) #sig.order有时候会有缺失???
group <- predict(estim.nmf2) # 提出亚型
table(group)
write.table(group,file="group4.txt",sep="\t",row.names=F,quote=F)

pdf("consensusmap4.pdf",onefile = FALSE,
    width = 8,
    height =8)
consensusmap(estim.nmf2,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(NMF.Exp.rank)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
dev.off()

pdf("basismap4.pdf",onefile = FALSE,
    width = 8,
    height =8)
basismap(estim.nmf2,
         cexCol = 1.5,
         cexRow = 1,
         annColors=list(c("deepskyblue3","lightgoldenrod2",'darkred','green4')))
dev.off()


library(tidyverse)
groupid=read.table("group4.txt",sep="\t",header=T,check.names=F)
expOSid=cbind(rt1,groupid=groupid$x)

library(survival)
library(survminer)
sfit <- survfit(Surv(OS.time, OS) ~ groupid,
                data = expOSid)

pdf("surv4.pdf",onefile = FALSE,
    width = 8,
    height =8)
ggsurvplot(sfit,pval = T,palette = c("deepskyblue3","lightgoldenrod2",'darkred','green1'))
dev.off()

library(Rtsne)
tsne_out = Rtsne(t(rt),perplexity = 30)
pdat = data.frame(tsne_out$Y,factor(group))
colnames(pdat) = c("Y1","Y2","group")
head(pdat)

library(ggplot2)
library(paletteer)
pdf("tsne4.pdf",onefile = FALSE,
    width = 8,
    height =8)
ggplot(pdat,aes(Y1,Y2))+
  geom_point(aes(Y1,Y2,fill = group),shape = 21,color = "black")+
  stat_ellipse(aes(color = group,fill = group),
               geom = "polygon",
               alpha = 0.3,
               linetype = 2)+
  scale_color_paletteer_d("RColorBrewer::Set3")+
  scale_fill_paletteer_d("RColorBrewer::Set3")+
  theme_classic()+
  theme(legend.position = "top")
dev.off()

groupid=read.table("group4.txt",sep="\t",header=T,check.names=F)
rt2=read.table("gene606OS1_log2.txt",sep="\t",header=T,check.names=F)
rtgroup4=data.frame(rt2,groupid)
colnames(rtgroup4)
ovca427group4=rtgroup4[,c(1,610)]
diffscal=read.table("diffmRNAExplog2Scale.txt",sep="\t",header=T,check.names=F,row.names = 1)
diffscalOvca=diffscal[,-c(1:88)]
diffscalOvca=as.data.frame(t(diffscalOvca))
diffscalOvca=cbind(sample=rownames(diffscalOvca),diffscalOvca)
diffOvcagroup4=inner_join(diffscalOvca,ovca427group4,by='sample')
diffOvcagroup4$x
diffOv4.mat=diffOvcagroup4
rownames(diffOv4.mat)=diffOv4.mat[,1]
diffOv4.mat=diffOv4.mat[,-c(1,ncol(diffOv4.mat))]
diffOv4.mat=t(diffOv4.mat)
#limma找差异基因
library(limma)
# 设定分组
condition <- factor(diffOvcagroup4$x)
condition
table(condition)
design1=model.matrix(~factor(condition))
fit1=lmFit(diffOv4.mat,design1)
fit1=eBayes(fit1)
options(digits = 4)
b2 <- topTable(fit1,coef=2,adjust='BH')
b3 <- topTable(fit1,coef=3,adjust='BH')
b4 <- topTable(fit1,coef=4,adjust='BH',n=Inf)
length(which(b4$adj.P.Val < 0.01))
top=b4[b4$adj.P.Val < 0.01,]
top3904list=rownames(top)
write.table(top3904list,file="top3904list.txt",sep="\t",quote=F,col.names=F,row.names = F) 
#1. raw count —edgeR—606genes—NMF 4 group
#新的流程，和ICGC的表达基因取交集后得到3583个基因后再建模

#下面的流程暂不用，仅保留做参考
#单因素回归分析
#install.packages('survival')
library(survival)
diffscal=read.table("diffmRNAExplog2Scale.txt",sep="\t",header=T,check.names=F,row.names = 1)
diffscalOvca=diffscal[,-c(1:88)]
diffscalOvca=as.data.frame(t(diffscalOvca))
diffscalOvca=cbind(sample=rownames(diffscalOvca),diffscalOvca)

OVCAid=read.table('OVCAidTCGA.txt',header = F,sep = '\t')
surdata=read.table('TCGA_survival_data.txt', header=T, sep="\t")
rownames(surdata)=surdata$sample
surdataOVCA=surdata[as.vector(OVCAid[,1]),]
surdataOVCAos=surdataOVCA[,c(1,3,2)]
table(which(is.na(surdataOVCAos[,3])))

mergeOS <- inner_join( surdataOVCAos,diffscalOvca,
                      by = "sample") 
mergeOS[1:4,1:4]

table(which(is.na(mergeOS)))
write.table(mergeOS,file="diffmRNAExplog2ScaleOS.txt",sep="\t",quote=F,col.names=T,row.names = F) 

gene3904=mergeOS[,colnames(mergeOS) %in% top3904list]
gene3904OS=data.frame(mergeOS[,1:3],gene3904)
rownames(gene3904OS)=gene3904OS[,1]
gene3904OS=gene3904OS[,-1]

outTab=data.frame()
for(i in colnames(gene3904OS[,3:ncol(gene3904OS)])){
  cox <- coxph(Surv(OS.time, OS) ~ gene3904OS[,i], data = gene3904OS)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)

sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<0.01,] #P鍊肩殑绛涢?夐槇鍊煎彲璁句负0.10
write.table(sigTab,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)

sigGenes=c("OS.time","OS")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=gene3904OS[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

#lasso回归
#install.packages("glmnet")
#install.packages("survival")
library("glmnet")

rt=read.table("uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)
which(is.na(rt))
rt=rt[-which(is.na(rt)),]
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$OS.time,rt$OS))

fit <- glmnet(x, y, family = "cox", maxit = 2000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 2000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("OS.time","OS",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp.txt",sep="\t",row.names=F,quote=F)

#多因素回归分析

rt=read.table("lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
multiCox=coxph(Surv(OS.time, OS) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)

#计算出每个患者的风险评分，并按中位数分为高低危
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("OS.time","OS",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk.txt",
            sep="\t",
            quote=F,
            row.names=F)
#K-M survival
#install.packages('survminer')
library(survival)
library(survminer)
rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(OS.time, OS) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(OS.time, OS) ~ risk, data = rt)
pdf(file="survival.pdf",onefile = FALSE,
    width = 8,
    height =8)
ggsurvplot(fit, 
           data=rt,
           #conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=T,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(days)",
           break.time.by = 365,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.18,
           risk.table.y.text = F,
           ncensor.plot = T,
           palette=c('indianred1','skyblue1'),
           ncensor.plot.height = 0.18,
           conf.int.style = "ribbon")
dev.off()
#summary(fit)

#RiskPlot
#install.packages("pheatmap")

rt=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10

#RiskScore
pdf(file="RiskScore.pdf",width = 12,height = 5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("skyblue1",lowLength),
           rep("indianred1",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
dev.off()

#SurvStat
color=as.vector(rt$OS)
color[color==1]="indianred1"
color[color==0]="skyblue2"
pdf(file="SurvStat.pdf",width = 12,height = 5)
plot(rt$OS.time,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (days)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

#heatmap  两种方法

rt1=rt[c(3:(ncol(rt)-1))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="hist.pdf",width = 12,height = 5)
hist(rt1)
dev.off()

#pheatmap 没法分组聚类
library(pheatmap)
bk <- c(seq(-2,2,by=0.2))
color = colorRampPalette(colors = c('blue','white',"red"))(100)
pdf(file="Hm-median.pdf",width = 15,height = 6)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = F,
         clustering_method = 'median',
         color = color,
         #breaks = bk,
         fontsize_row=8,
         fontsize_col=3,
         show_colnames=F)
dev.off()

#分组聚类！！！
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

table(annotation$type)
Group=factor(annotation$type)
col_fun = colorRamp2(c(-1.5, 0, 1.5),c("#2fa1dd", "white", "red"))
#col = colorRamp2(breaks = seq(min(rt1), max(rt1), length = 3), colors = c("#2fa1dd", "white", "#f87669"))

top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c( "#2fa1dd","#f87669")),
                       labels = c("low","high"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

pdf(file="new-Hm1.5.pdf",width = 15,height = 6)
Heatmap(rt1,name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = Group,
        show_heatmap_legend = T,
        border = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL)
dev.off()


##ROC
#install.packages("survivalROC")

library(survivalROC)
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)

#1年ROC
pdf(file="ROC-1.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =365, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#3年ROC
pdf(file="ROC-3.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =365*3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#5年ROC
pdf(file="ROC-5.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =365*5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#10年ROC
pdf(file="ROC-10.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =365*10, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#14年ROC
pdf(file="ROC-14.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =365*14, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#整合1，3，5,10,14年ROC
pdf(file="ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)

roc1=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                 predict.time =365, method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("riskScore ROC curve"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

roc2=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                 predict.time =365*3, method="KM")   #在此更改时间，单位为年
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="maroon1",lwd=2)
#text(locator(1), paste("1 year",round(roc2$AUC,3),sep=":"),col="blue")

roc3=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                 predict.time =365*5, method="KM")   #在此更改时间，单位为年
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="yellow4",lwd=2)
#text(locator(1), paste("2 year",round(roc2$AUC,3),sep=":"),col="green")

roc4=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                 predict.time =365*10, method="KM")   #在此更改时间，单位为年
lines(roc4$FP,roc4$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="olivedrab3",lwd=2)
#text(locator(1), paste("2 year",round(roc2$AUC,3),sep=":"),col="yellow")

roc5=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                 predict.time =365*14, method="KM")   #在此更改时间，单位为年
lines(roc5$FP,roc5$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="tomato",lwd=2)
#text(locator(1), paste("1 year",round(roc2$AUC,3),sep=":"),col="blue")

legend("bottomright", 
       c("1-year AUC:0.780","3-year AUC:0.748","5-year AUC:0.706","10-year AUC:0.908","14-year AUC:0.968"),
       lwd=2,
       col=c("red","maroon1","yellow4","olivedrab3","tomato"))
dev.off()

##C-index
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("survcomp")

library(survcomp)
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
cindex <- concordance.index(x=rt$riskScore,
                            surv.time = rt$OS.time, 
                            surv.event = rt$OS,
                            method = "noether")
cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

library(LPS)
#LPS predict system

rt=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])

rt1=rt[,c(3:(ncol(rt)-2))]
rtexp=t(t(rt1))

rt.cli=rt[c(1,2,ncol(rt))]
group <- rt.cli$risk 
heat.map(rtexp)
coeff <- LPS.coeff(data=rtexp, response=group) #from t test value
m <- LPS(data=rtexp, coeff=coeff, response=group, k=23)
plot(m, "probability", yaxt="s")
predict(m, rtexp, type="probability", plot=TRUE)

# from multicox coeff
outTabmulti=as.data.frame(outTab)
multicoeff=outTabmulti[,c(1,2,6)]
rownames(multicoeff)=multicoeff[,1]
multicoeff=multicoeff[,-1]
colnames(multicoeff)[1]='t'
dimnames=list(rownames(multicoeff),colnames(multicoeff))
multicoeff=matrix(as.numeric(as.matrix(multicoeff)),nrow=nrow(multicoeff),dimnames=dimnames)

m <- LPS(data=rtexp, coeff=multicoeff, response=group, k=23)
plot(m, "probability", yaxt="s")
pdf('multicoeffLPS.pdf',width = 10,height = 12)
predict(m, rtexp, type="probability",threshold = 0.8, plot=TRUE)
dev.off()

predict(m, rtexp, type="score",threshold = 0.8, plot=TRUE)
# Class prediction plot 
predict(m, rtexp, plot=TRUE) 
# Wright et al. class prediction 
table( group, prediction = predict(m, rtexp), exclude = NULL ) 
# More stringent threshold 
table( group, prediction = predict(m, rtexp, threshold=0.99), exclude = NULL ) 
# Radmacher et al. class prediction 
table( group, prediction = predict(m, rtexp, method="Radmacher"), exclude = NULL ) 
# Probabilities 
predict(m, rtexp, type="probability", method="Wright") 
predict(m, rtexp, type="probability", method="Radmacher") 
predict(m, rtexp, type="probability", method="exact") 
# Probability plot 
predict(m, rtexp, type="probability", plot=TRUE) 
# Annotated probability plot 
side <- data.frame(group, row.names=rownames(rtexp)) 
predict(m, rtexp, side=side, type="probability", plot=TRUE) 
# Score plot 
predict(m, rtexp, type="score", plot=TRUE)

#install.packages("timeROC")


#引用包
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

