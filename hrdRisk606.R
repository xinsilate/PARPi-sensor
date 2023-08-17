#if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
#BiocManager::install("rtracklayer")

rm(list=ls())
setwd("/Users/sunx/Documents/ovca-immune/ucsc/NMF606/HRD/")

gz=gzfile('/Users/sunx/Documents/ovca-immune/ucsc/explore/HRD/TCGA.HRD_withSampleID.txt.gz','rt')
data=read.table(gz,sep = "\t",header = T,check.names = F)
data[1:4,1:6]

names=colnames(data)
types=read.table('TCGA_GTEX_category.txt',sep = '\t',header = T,check.names = F)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
data[1:6,1:6]
#ovary=types %>% filter(str_detect(TCGA_GTEX_main_category,'Ovary'))
ovarian=types %>% filter(str_detect(TCGA_GTEX_main_category,'Ovarian'))
#ovaryID=c(ovary[,1])
#ovarydata=select(data,ovaryID)
ovarianID=c('sampleID',ovarian[,1])
interID=intersect(names,ovarianID)
hrddata=select(data,interID)
rownames(hrddata)=hrddata[,1]
hrddata=hrddata[,-1]
hrddata=t(hrddata)
hrddata=as.data.frame(hrddata)
hrddata=cbind(id=rownames(hrddata),hrddata)

riskdata=read.table('risk23.txt',sep = '\t',header = T,check.names = F)

hrdrisk <- inner_join(hrddata, riskdata,
                      by = "id") 
as.factor(hrdrisk$risk)
write.table(hrdrisk, file="hrdrisk.txt", sep="\t", quote=F, row.names=F, col.names=T)

#重新读入文件
hrdrisk=read.table('hrdrisk.txt',sep = '\t',header = T,check.names = F,row.names = 1)

#把数据转换成ggplot2输入文件
data1=hrdrisk
data1=data1[,c(1:4,31)]
data=melt(data1, id.vars=c("risk"))
colnames(data)=c("Type", "Gene", "Expression")
data$Type=factor(data$Type,levels = c('low','high'))

#绘制箱线图

p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            ylab="Gene expression",
            xlab="",
            add='jitter',
            legend.title="Type",
            palette = c("skyblue2", "indianred2"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出箱线图
pdf(file="boxplot.hrd.pdf", width=7, height=5)
print(p1)
dev.off()


pdf(file = "HRD.pdf")
# Plot the chart.
boxplot(HRD ~ risk, data = hrdrisk, xlab = "risk",
        ylab = "HRD", main = "HRD")
dev.off()
wilcox.test(HRD ~ risk,data = hrdrisk)


pdf(file = "hrd-loh.pdf")
# Plot the chart.
boxplot(`hrd-loh` ~ risk, data = hrdrisk, xlab = "risk",
        ylab = "hrd-loh", main = "hrd-loh")
dev.off()
wilcox.test(`hrd-loh` ~ risk,data = hrdrisk)

pdf(file = "ai1.pdf")
# Plot the chart.
boxplot(ai1 ~ risk, data = hrdrisk, xlab = "risk",
        ylab = "ai1", main = "ai1")
dev.off()
wilcox.test(ai1 ~ risk,data = hrdrisk)


pdf(file = "lst1.pdf")
# Plot the chart.
boxplot(lst1 ~ risk, data = hrdrisk, xlab = "risk",
        ylab = "lst1", main = "lst1")
dev.off()
wilcox.test(lst1 ~ risk,data = hrdrisk)


#重新读入文件
hrdrisk=read.table('hrdrisk.txt',sep = '\t',header = T,check.names = F,row.names = 1)
data=hrdrisk
#获取肿瘤突变负荷最优的cutoff
res.cut=surv_cutpoint(data, time = "OS.time", event = "OS", variables =c("HRD"))
cutoff=as.numeric(res.cut$cutpoint[1])
hrdType=ifelse(data[,"HRD"]<=cutoff, "L-HRD", "H-HRD")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(hrdType, "+", scoreType)

#定义生存分析函数,绘制肿瘤突变负荷的生存曲线
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(OS.time, OS) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time, OS) ~ group, data = surData)
  #print(surv_median(fit))
  
  #绘制生存曲线
  bioCol=c("skyblue4","indianred4","skyblue1","indianred1","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     #surv.median.line = "hv",
                     risk.table=F,
                     cumevents=F,
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, width=5.5, height=4.8, onefile = FALSE)
  print(surPlot)
  dev.off()
}

#调用函数,绘制肿瘤突变负荷的生存曲线
data$group=hrdType
bioSurvival(surData=data, outFile="HRD.survival.pdf")


#定义生存分析函数,绘制肿瘤突变负荷联合病人风险的生存曲线
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(OS.time, OS) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time, OS) ~ group, data = surData)
  #print(surv_median(fit))
  
  #绘制生存曲线
  bioCol=c("indianred4","skyblue4","indianred1","skyblue1","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     #surv.median.line = "hv",
                     risk.table=F,
                     cumevents=F,
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, width=5.5, height=4.8, onefile = FALSE)
  print(surPlot)
  dev.off()
}

#调用函数,绘制肿瘤突变负荷联合病人风险的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="HRD-risk.survival.pdf")



