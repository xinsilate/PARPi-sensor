
rm(list=ls())
setwd("/Users/sunx/Documents/ovca-immune/ICGC/")

library(data.table)
library(dplyr)
library(tibble)
library(tidyr)

###exp表达文件的提取和长宽转换,提取raw read count,不清楚标准化方法
expr <- fread('exp_seq.OV-AU.tsv.gz',data.table = F)%>%
  dplyr::select(icgc_donor_id,gene_id,raw_read_count)%>%
  group_by(icgc_donor_id,gene_id)%>%
  summarise_all(max)%>%
  pivot_wider(names_from = 'icgc_donor_id',values_from = "raw_read_count")%>%
  summarise_all(function(x){ifelse(is.na(x),0,x)})

#normalize数据
library(limma)
expr[1:4,1:4]

#这个没有重复的ensemble ID，去重没有意义
rt=distinct(expr,gene_id,.keep_all = T)
write.table(rt, file="OVAU.ensem.rawCount.txt", row.names=F, col.names=T, quote=F,sep="\t")
rt=as.matrix(expr)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#rt=avereps(rt)
rt1=normalizeBetweenArrays(as.matrix(rt))
rta=cbind(id=rownames(rt1),rt1)
write.table(rta, file="OVAU.ensem.limmaNorm.txt", row.names=F, col.names=T, quote=F,sep="\t")


#方法1，重新读入数据
OVAU.mat0=read.table('OVAU.ensem.limmaNorm.txt',header=T,sep='\t',check.names = F)
#if need, remove number behind the dot of Ensembl_ID
#COADdata1<-separate(COADdata,Ensembl_ID,into= c("Ensembl_ID"),sep="\\.")
OVAU.mat.data=OVAU.mat0
transfer=read.table('probeMap_gencode.v23.annotation.gene.probemap.txt',header=T,sep='\t',check.names = F)
# remove number behind the dot of Ensembl_ID
transfer<-separate(transfer,id,into= c("id"),sep="\\.")
which(transfer$id %in% 'ENSG00000275791')
OVAU.symb3=inner_join(transfer,OVAU.mat.data,by='id')

#去重的时候只保留最大的那个其实保持了数据的随机性，适用性更强
OVAU.symb4=distinct(OVAU.symb3,gene,.keep_all = T)
colnames(OVAU.symb4)[1:8]
OVAU.symb4=OVAU.symb4[,-c(1,3,4,5,6)]
#删除gene为NA的行
which(is.na(OVAU.symb4$gene))

rownames(OVAU.symb4)=OVAU.symb4[,1]
OVAU.symb4=OVAU.symb4[,-1]
OVAU.symb4=t(t(OVAU.symb4))
hist(OVAU.symb4)
which(rowMeans(OVAU.symb4)=='0')
OVAU.symb4=OVAU.symb4[rowMeans(OVAU.symb4)>0,]
OVAU.symb4=log2(OVAU.symb4+1)
hist(OVAU.symb4)
OVAU.symb4a=cbind(id=rownames(OVAU.symb4),OVAU.symb4)
write.table(OVAU.symb4a, file="OVAU.symb.limmaNormlog2.txt", row.names=F, col.names=T, quote=F,sep="\t")
dev.off()

gene3904=read.table('top3904list.txt',sep = '\t',check.names = F,header = F)
inter=intersect(gene3904$V1,rownames(OVAU.symb4))
write.table(inter, file="top3583list.txt", row.names=F, col.names=F, quote=F,sep="\t")
which(gene3904 %in% 'RP11.678G14.3')

#新的流程 cox model : raw count— limma norm—log2—scale—model
#单因素回归分析
#install.packages('survival')
library(survival)

ovcaRaw0=read.table('ovcaGTEx516GeneCountlog2.txt',header=T,sep='\t',check.names = F,row.names = 1)
#去正常卵巢组织
ovcaRaw=ovcaRaw0[,-c(1:88)]
#去log化,edgeR DEG
ovcaRaw=2^ovcaRaw-1
#ovcaRaw=ovcaRaw[rowMeans(ovcaRaw)>1,]
#limma norm
ovcaRaw=normalizeBetweenArrays(as.matrix(ovcaRaw))
#log2
ovcaRaw=log2(ovcaRaw+1)
ovcaRawlog2a=ovcaRaw
ovcaRawlog2=cbind(id=rownames(ovcaRaw),ovcaRaw)
write.table(ovcaRawlog2, file="ovca427.raw.limmaNorm.log2.txt", row.names=F, col.names=T, quote=F,sep="\t")

#scale 去中心化
library(vegan) #标准化和中心化
ovcaRaw=decostand(ovcaRaw,'standardize',MARGIN = 1)
ovcaRawa=cbind(id=rownames(ovcaRaw),ovcaRaw)
write.table(ovcaRawa, file="ovca427.raw.limmaNorm.log2.scale.txt", row.names=F, col.names=T, quote=F,sep="\t")

OVCAid=read.table('OVCAidTCGA.txt',header = F,sep = '\t')
surdata=read.table('TCGA_survival_data.txt', header=T, sep="\t")
rownames(surdata)=surdata$sample
surdataOVCA=surdata[as.vector(OVCAid[,1]),]
surdataOVCAos=surdataOVCA[,c(1,3,2)]
surdataOVCAos=surdataOVCAos[-which(is.na(surdataOVCAos[,2])),]
table(which(is.na(surdataOVCAos[,2])))
surdataOVCAos$OS.time=surdataOVCAos$OS.time/365
colnames(surdataOVCAos)[1]='id'
write.table(surdataOVCAos, file="ovca425.OS.txt", row.names=F, col.names=T, quote=F,sep="\t")

ovcaRaw=t(ovcaRaw)
ovcaRawb=cbind(id=rownames(ovcaRaw),ovcaRaw)
ovcaRawb=as.data.frame(ovcaRawb)
mergeOS <- inner_join(surdataOVCAos, ovcaRawb,
                      by = "id") 
write.table(mergeOS, file="ovca425.OS.raw.limmaNorm.log2.scale.txt", row.names=F, col.names=T, quote=F,sep="\t")

ovcaRawlog2a=t(ovcaRawlog2a)
ovcaRawlog2a=cbind(id=rownames(ovcaRawlog2a),ovcaRawlog2a)
ovcaRawlog2a=as.data.frame(ovcaRawlog2a)
mergeOSa <- inner_join(surdataOVCAos, ovcaRawlog2a,
                      by = "id") 
write.table(mergeOSa, file="ovca425.OS.raw.limmaNorm.log2.txt", row.names=F, col.names=T, quote=F,sep="\t")
gene606=read.table('gene606.txt',header = F,sep = '\t')
gene606=c('id','OS.time','OS',gene606$V1)
gene606os.mat=mergeOSa[,colnames(mergeOSa) %in% gene606]
write.table(gene606os.mat, file="ovca425.606geneOS.raw.limmaNorm.log2.txt", row.names=F, col.names=T, quote=F,sep="\t")

#重新读入数据
ovca425OS=read.table('ovca425.OS.raw.limmaNorm.log2.scale.txt',sep = '\t',check.names = F,header = T,row.names = 1)
gene3583list=read.table('top3583list.txt',sep = '\t',check.names = F,header = F)
gene3583list=c('OS.time','OS',gene3583list$V1)
gene3583os=ovca425OS[,colnames(ovca425OS) %in% gene3583list]
head(colnames(gene3583os))


outTab=data.frame()
for(i in colnames(gene3583os[,3:ncol(gene3583os)])){
  cox <- coxph(Surv(OS.time, OS) ~ gene3583os[,i], data = gene3583os)
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
write.table(outTab,file="3583uniCoxResult.limma.txt",sep="\t",row.names=F,quote=F)

sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<0.01,] #P鍊肩殑绛涢?夐槇鍊煎彲璁句负0.10
write.table(sigTab,file="3583uniCoxResult.Sig169.limma.txt",sep="\t",row.names=F,quote=F)

sigGenes=c("OS.time","OS")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=gene3583os[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="3583uniSigExp169.limma.txt",sep="\t",row.names=F,quote=F)

#lasso回归
#install.packages("glmnet")
#install.packages("survival")
library("glmnet")

rt=read.table("3583uniSigExp169.limma.txt",header=T,sep="\t",row.names=1,check.names=F)
which(is.na(rt))
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$OS.time,rt$OS))

fit <- glmnet(x, y, family = "cox", maxit = 5000)
pdf("lambda169.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 5000)
pdf("cvfit169.pdf")
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
write.table(lassoSigExp,file="169lassoSigExp45.txt",sep="\t",row.names=F,quote=F)

#多因素回归分析

rt=read.table("169lassoSigExp45.txt",header=T,sep="\t",check.names=F,row.names=1)
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
write.table(outTab,file="45multiCox23.txt",sep="\t",row.names=F,quote=F)

#计算出每个患者的风险评分，并按中位数分为高低危
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("OS.time","OS",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk23.txt",
            sep="\t",
            quote=F,
            row.names=F)
#K-M survival
library(survival)
library(survminer)
rt=read.table("risk23.txt",header=T,sep="\t")
diff=survdiff(Surv(OS.time, OS) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(OS.time, OS) ~ risk, data = rt)
pdf(file="survival.23.pdf",onefile = FALSE,
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
           xlab="Time(years)",
           break.time.by = 1,
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















#和生存数据合并
sur.OVAU.org <- read_delim("donor.OV-AU.tsv.gz", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
colnames(sur.OVAU.org)
sur.OVAU=sur.OVAU.org[,c('icgc_donor_id','donor_survival_time','donor_vital_status','donor_age_at_diagnosis','donor_tumour_stage_at_diagnosis')]
write.table(sur.OVAU, file="sur.OVAU.txt", row.names=F, col.names=T, quote=F,sep="\t")
colnames(sur.OVAU)=c('id','OS.time','OS','age','stage')
sur.OVAU$OS[sur.OVAU$OS=='deceased']=0
sur.OVAU$OS[sur.OVAU$OS=='alive']=1

OVAU.symb.org=read.table('OVAU.symb.limmaNormlog2.txt',header=T,sep='\t',check.names = F,row.names = 1)
library(vegan) #标准化和中心化
OVAU.symb=decostand(OVAU.symb.org,'standardize',MARGIN = 1)
OvauSymb=t(OVAU.symb.org)
OvauSymb=as.data.frame(OvauSymb)
OvauSymb=cbind(id=rownames(OvauSymb),OvauSymb)
which(colnames(OvauSymb) %in% 'RP11.678G14.3')

write.table(OvauSymb, file="OVAU.symb.limmaNormlog2.scale.txt", row.names=F, col.names=T, quote=F,sep="\t")

OvauOS=inner_join(sur.OVAU,OvauSymb,by='id')
write.table(OvauOS, file="OvauOS.limmaNormlog2.scale.txt", row.names=F, col.names=T, quote=F,sep="\t")

#对新数据计算riskScore
#计算出每个患者的风险评分，并按中位数分为高低危
### check.names=T的作用是把基因中的横线-变为小点.。
OvauOS.data=read.table("OvauOS.limmaNormlog2.scale.txt",header=T,sep="\t",check.names=F,row.names=1)
riskScore=predict(multiCox,type="risk",newdata=OvauOS.data)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("OS.time","OS",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(OvauOS.data[,outCol],riskScore,risk)),cbind(OvauOS.data[,outCol],riskScore,risk)),
            file="risk.ovau.txt",
            sep="\t",
            quote=F,
            row.names=F)
#K-M kurve for OV-AU
risk.ovau=read.table("risk.ovau.txt",header=T,sep="\t",check.names=F,row.names=1)
rt=risk.ovau
diff=survdiff(Surv(OS.time, OS) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(OS.time, OS) ~ risk, data = rt)
pdf(file="survival23.ovauPre.pdf",onefile = FALSE,
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

