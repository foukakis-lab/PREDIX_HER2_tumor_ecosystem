#Fig.a1
library(ggpubr);library(data.table)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
table(clin$Arm)
clin$pCR[clin$Arm=="DHP"]%>%table()
clin$pCR[clin$Arm=="T-DM1"]%>%table()
# Data
df <- data.frame(Arm=c("DHP", "T-DM1"),
                 pCR=c(45.5,43.9))
ggbarplot(df, "Arm", "pCR",
          fill = "Arm", color = "white",
          palette = c("#8491B4FF","#91D1C2FF"),
          label = TRUE, lab.pos = "in", lab.col = "white")
# Portrait 3X5
median(clin$OS.time)
library("survival");library(survminer)
sfit <- survfit(Surv(EFS.time,EFS.status)~Arm, data=clin)
sfit
summary(sfit)

fit <- survfit(Surv(EFS.time,EFS.status) ~ Arm, data = clin)
ggsurvplot(fit, palette = c("#8491B4FF","#91D1C2FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(Arm)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES)+as.factor(pCR), data=clin) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))

#Fig.a2
library(ComplexUpset);library(ggvenn);library(ggsci);library("scales");library(ggplot2);library(readxl);library(data.table)
show_col(pal_npg("nrc")(10))
show_col(pal_npg("nrc", alpha = 0.6)(10))
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
meta=meta[,c("patientID","Arm","ER")]
# merge
meta$patientID=as.character(meta$patientID)
wes$patientID=as.character(wes$patientID)
swgs$patientID=as.character(swgs$patientID)
rna$patientID=as.character(rna$patientID)
MS$patientID=as.character(MS$patientID)
images$patientID=as.character(images$patientID)
data=left_join(meta,wes,by="patientID")%>%left_join(swgs,by="patientID")%>%
  left_join(rna,by="patientID")%>%left_join(MS,by="patientID")%>%
  left_join(images,by="patientID")%>%left_join(xenium,by="patientID")
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
table(data$xenium)
write.csv(data, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Fig1a.csv")
id=c("WES","sWGS","RNAseq","MS","Image","xenium")
library(ComplexHeatmap)
library(circlize)
data=data[order(data$Arm),]
ha = HeatmapAnnotation(Arm=data$Arm,
                       ER = data$ER,
                       WES = data$WES,
                       sWGS=data$sWGS,
                       RNAseq=data$RNAseq,
                       MS=data$MS,
                       Image=data$Image,
                       Xenium=data$xenium,
                       col = list(Arm =c("DHP" =  "#8491B4FF", "T-DM1" = "#91D1C2FF"),
                                  ER=c("positive" =  "#616569", "negative" = "#eeeeee"),
                                  WES=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  sWGS=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  RNAseq=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  DNA_cycle2=c("1" = "#916ba6", "0" = "#eeeeee"),
                                  MS=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  Image=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  Xenium=c("1" =  "#916ba6", "0" = "#eeeeee")))


draw(ha)

set.seed(123)
mat = matrix(rnorm(80, 2), 8, 197)
Heatmap(mat, cluster_columns=FALSE,top_annotation = ha)
# 11X6 P