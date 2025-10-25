############Upset plot###############
library(readxl);library(data.table);library(tidyverse)
wes=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.txt")
wes=data.frame(sampleID.wes=unique(wes$Tumor_Sample_Barcode),WES=1)%>%mutate(patientID=substr(sampleID.wes,9,12))
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
rna=data.frame(sampleID.rna=colnames(rna),RNAseq=1)%>%mutate(patientID=substr(sampleID.rna,9,12))
swgs=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_amp_peak_curated.txt")
swgs=data.frame(patientID=swgs$patientID,sWGS=1)
images=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/image_metrics_PREDIX_HER2.rds")
images$patientID=as.character(images$patientID)
images$Image=1
images=images[,c("patientID","Image")]
# MS
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
MS$sampleID.ms=MS$Code
MS$MS=1
MS=MS[,c("sampleID.ms","patientID","MS")]
MS=MS[!duplicated(MS$patientID),]
# Xenium
xenium=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/xenium_meta.csv")
xenium=xenium[xenium$tpt=="pre",]
xenium$xenium=1
xenium=xenium[,c("patientID","xenium")]
xenium$patientID=as.character(xenium$patientID)
# clin
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta=meta[,c("patientID","Arm","ER","ANYNODES","TUMSIZE","Response")]
# merge
meta$patientID=as.character(meta$patientID)
wes$patientID=as.character(wes$patientID)
swgs$patientID=as.character(swgs$patientID)
rna$patientID=as.character(rna$patientID)
MS$patientID=as.character(MS$patientID)
images$patientID=as.character(images$patientID)
data=left_join(meta,wes,by="patientID")%>%left_join(swgs,by="patientID")%>%left_join(rna,by="patientID")%>%
  left_join(MS,by="patientID")%>%left_join(images,by="patientID")%>%left_join(xenium,by="patientID")
data$WES[is.na(data$WES)]=0
data$sWGS[is.na(data$sWGS)]=0
data$RNAseq[is.na(data$RNAseq)]=0
data$MS[is.na(data$MS)]=0
data$Image[is.na(data$Image)]=0
data$xenium[is.na(data$xenium)]=0
table(data$WES)
table(data$sWGS)
table(data$RNAseq)
table(data$MS)
table(data$Image)
data=as.data.frame(data)
id=c("WES","sWGS","RNAseq","MS","Image","xenium")
library(UpSetR)
upset(data, sets = id, sets.bar.color = "#56B4E9",
      order.by = "freq")
# 5X5 
##########Export meta data########
# Read data
data=data[,c("patientID","sampleID.wes","sampleID.rna","Arm","ER","ANYNODES","TUMSIZE","Response")]
# CUTseq
meta = as.data.frame(fread("E:/Projects/PREDIX_HER2/CUTseq/data/meta/meta_CUTseq.txt")) 
meta=meta[meta$Sample_Cohort=="PREDIX HER2"&meta$tpt=="Baseline",]
cutseq=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_gene_curated.txt")%>%as.data.frame()
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count$patientID=substr(seg_count$sample,1,4)
pid=intersect(seg_count$patientID,cutseq$patientID)
seg_count=seg_count[seg_count$patientID%in%pid,]
meta=meta[meta$sampleID%in%seg_count$sample,]
data[is.na(data)]=NA
meta=meta[,c("sampleID","patientID","BARCODE")]
colnames(meta)=c("sampleID.cutseq","patientID","barcode_cutseq")
meta$patientID=as.character(meta$patientID)
data=left_join(data,meta,by="patientID")
data=data[,c("patientID","sampleID.wes","sampleID.rna","sampleID.cutseq","barcode_cutseq","Arm","ER","ANYNODES","TUMSIZE","Response")]
write.table(data,file="E:/Projects/PREDIX_HER2/Multimodal/Data_repository/PREDIX_HER2_multiomics_meta.txt",quote = F,row.names =F,sep="\t")
########### FigureS1.B ###########
library("survival");library(survminer)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
fit <- survfit(Surv(EFS.time,EFS.status) ~ Response, data = data)
summary(fit)

fit <- survfit(Surv(EFS.time,EFS.status) ~ Response, data = data)
ggsurvplot(fit, data =data, palette = c("#fdb462","#20A39E"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(data$Response)
# 4X4
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(Response)+Arm+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=data) ##DII_density_with_supp
ShowRegTable(cox.test)

DHP=data[data$Arm=="DHP",]
fit <- survfit(Surv(EFS.time,EFS.status) ~ Response, data = DHP)
summary(fit)
ggsurvplot(fit, data =DHP, palette = c("#fdb462","#20A39E"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(DHP$Response)
# 4X4
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(Response)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=DHP) ##DII_density_with_supp
ShowRegTable(cox.test)

TDM1=data[data$Arm=="T-DM1",]
fit <- survfit(Surv(EFS.time,EFS.status) ~ Response, data = TDM1)
summary(fit)
ggsurvplot(fit, data =TDM1, palette = c("#fdb462","#20A39E"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(TDM1$Response)
# 4X4
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(Response)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=TDM1) ##DII_density_with_supp
ShowRegTable(cox.test)





library("survival");library(survminer)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
fit <- survfit(Surv(EFS.time,EFS.status) ~ Arm, data = data[data$HER2neu1=="2+",])
summary(fit)

fit <- survfit(Surv(EFS.time,EFS.status) ~ Arm, data = data[data$HER2neu1=="2+",])
ggsurvplot(fit, data =data[data$HER2neu1=="2+",], palette = c("#fdb462","#20A39E"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(data$Response)

fit <- survfit(Surv(EFS.time,EFS.status) ~ Arm, data = data[data$HER2neu1=="3+",])
ggsurvplot(fit, data =data[data$HER2neu1=="3+",], palette = c("#fdb462","#20A39E"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
