################################
###########FigureS3b############
################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_gene_curated.txt")%>%as.data.frame()
drivers=fread("E:/Projects/Collaboration/BEVPAC/CUTseq/genelist_nik-zainal-etal.tsv") 
drivers=drivers$Gene
cna=cna[,c(drivers,"patientID")]
amp_gene=c("FGFR1","CCND1","MDM2","ERBB2","ZNF217","PIK3CA")
del_gene=c("MAP3K1","CDKN2A","CDKN2B","BRCA2","RB1","AKT1","MAP2K4","TP53","NCOR1")
cna[,del_gene]=-1*cna[,del_gene]
cna=cna[,c(del_gene,amp_gene,"patientID")]

genomic=left_join(cna,clin,by="patientID")%>%as.data.frame()
str(genomic)
colnames(genomic)
results=Logistic_batch_adjER(genomic,"pCR","Arm",c(del_gene,amp_gene),"ER")%>%as.data.frame()
results$group[results$biomarker%in%amp_gene]="amp"
results$group[results$biomarker%in%del_gene]="del"
data=results[,c("biomarker","group","whole_OR","Whole_LCI","Whole_UCI","whole_lr_p")]
data$whole_OR=as.numeric(data$whole_OR)
data$whole_lr_p=as.numeric(data$whole_lr_p)
data$FDR=p.adjust(data$whole_lr_p, method = "BH")
data$lnOR=data$whole_OR%>%log(base=exp(1))
data$log10FDR=-log10(data$FDR)

library(ggplot2)
library(ggrepel)
# FDR and color
thr_5FDR <- -log10(0.05)
color_map <- c("amp" = "brown", "del" = "steelblue")
# 绘图
p1 <- ggplot(data, aes(x = lnOR, y = log10FDR, fill= group)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.3,
             position = position_jitter(width = 0.05, height = 0.02)) +  # 外框点
  scale_fill_manual(values = color_map) +
  geom_hline(yintercept = thr_5FDR, color = "brown", linetype = "solid") +
  geom_vline(xintercept = 0, color = "gray40", linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(data, FDR < 0.25),
                           aes(label = biomarker),
                           size = 4, color = "black") +
  annotate("text", x = 2, y = thr_5FDR + 0.05, label = "5% FDR", color = "brown", hjust = 1) +
  annotate("text", x = -1, y = 0.2, label = "RD associated", color = "black", hjust = 0, size = 4) +
  annotate("text", x = 2, y = 0.2, label = "pCR associated", color = "black", hjust = 1, size = 4) +
  
  labs(x = "ln(ORadjusted)", y = expression(-log[10](FDR))) +
  xlim(-1, 2.5) +
  ylim(0, 2) +
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank())
p1

# 4X5
################################
###########FigureS3c############
################################
library(tidyverse);library(survminer);library(data.table);library(survival)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")%>%as.data.frame()
CNA=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/CUTseq_gene_baseline_complemental.rds")
mutation=genomic[,c("patientID","coding_mutation_TP53_oncokb")]
CNA=CNA[,c("BRCA2","FGFR1","NCOR1","ERBB2","PIK3CA","MAP3K1")]
colnames(CNA)=paste0("CNA_",colnames(CNA))
CNA$patientID=row.names(CNA)%>%as.integer()
data=left_join(clin,mutation,by="patientID")%>%left_join(CNA,by="patientID")%>%as.data.frame()
data$CNA_BRCA2_Del[data$CNA_BRCA2%in%c(0,1,2)]="No"
data$CNA_BRCA2_Del[data$CNA_BRCA2%in%c(-2,-1)]="Yes"
data$CNA_ERBB2_Amp[data$CNA_ERBB2%in%c(0,-1,-2)]="No"
data$CNA_ERBB2_Amp[data$CNA_ERBB2%in%c(2,1)]="Yes"
dhp=data[data$Arm=="DHP",]
tdm1=data[data$Arm=="T-DM1",]

#dhp
table(dhp$CNA_ERBB2_Amp)
fit <- survfit(Surv(EFS.time,EFS.status) ~ CNA_ERBB2_Amp, data = dhp)
ggsurvplot(fit,data = dhp, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(CNA_ERBB2_Amp)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=dhp) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))

# 4X4 
#tdm1 d,f,i
table(tdm1$CNA_ERBB2_Amp)
fit <- survfit(Surv(EFS.time,EFS.status) ~ CNA_ERBB2_Amp, data = tdm1)
ggsurvplot(fit,data = tdm1, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(CNA_ERBB2_Amp)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=tdm1) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))

table(tdm1$CNA_BRCA2_Del)
fit <- survfit(Surv(EFS.time,EFS.status) ~ CNA_BRCA2_Del, data = tdm1)
ggsurvplot(fit,data = tdm1, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(CNA_BRCA2_Del)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=tdm1) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))

#MSKCC -DHP e
library("survival");library(survminer);library(tableone);library(data.table);library(tidyverse)
tumor=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_clinical_sample.txt")
#tumor=tumor[tumor$SAMPLE_TYPE=="Primary",]
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_clinical_patient.txt")
meta=left_join(meta,tumor,by="PATIENT_ID")
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_cna.txt")%>%as.data.frame()
row.names(cna)=cna$Hugo_Symbol;cna$Hugo_Symbol=NULL;cna=cna["ERBB2",]
cna=t(cna)%>%as.data.frame()
cna$SAMPLE_ID=row.names(cna)
meta=left_join(meta,cna,by="SAMPLE_ID")
meta$PFS.STATUS=substr(meta$PFS_STATUS,1,1)
meta$PFS.STATUS=as.numeric(meta$PFS.STATUS)
meta=meta%>%filter(!is.na(meta$PFS.STATUS),PFS_MONTHS>3) 


fit <- survfit(Surv(PFS_MONTHS,PFS.STATUS) ~ ERBB2, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

table(meta$ERBB2)
cox.test <- coxph(Surv(PFS_MONTHS,PFS.STATUS)~as.factor(ERBB2)+as.factor(DX_PSTAGE)+SAMPLE_TYPE, data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)

# MSK ERBB2 CNA T-DM1 and T-dXD g,h
library("survival");library(survminer);library(tableone);library(data.table);library(tidyverse)
pts=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_clinical_patient.txt")
sample=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_clinical_sample.txt")
meta=left_join(sample,pts,by="PATIENT_ID")%>%filter(CANCER_TYPE=="Breast Cancer",HER2=="Yes",SAMPLE_TYPE=="Metastasis")%>%as.data.frame()
meta$OS_STATUS[meta$OS_STATUS=="0:LIVING"]=0
meta$OS_STATUS[meta$OS_STATUS=="1:DECEASED"]=1
meta$OS_STATUS=as.numeric(meta$OS_STATUS)
agent=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_timeline_treatment.txt")
agent$duration=agent$STOP_DATE-agent$START_DATE
#msk_chord_2024,T-MD1
adc_pts=agent$PATIENT_ID[agent$AGENT%in%c("ADO-TRASTUZUMAB EMTANSINE")&agent$START_DATE>0]
adc_pts=intersect(adc_pts,agent$PATIENT_ID[agent$AGENT!= c("PERTUZUMAB")])%>%unique()
adc_pts=intersect(meta$PATIENT_ID,adc_pts)
#tdm1_pts=intersect(agent$PATIENT_ID[agent$AGENT%in%c("TRASTUZUMAB")&agent$START_DATE<0&agent$duration>200],agent$PATIENT_ID[agent$AGENT%in%c("PERTUZUMAB")&agent$START_DATE<0&agent$duration>200])
meta=meta[meta$PATIENT_ID%in%adc_pts,]
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_cna.txt")%>%as.data.frame()
row.names(cna)=cna$Hugo_Symbol
cna$Hugo_Symbol=NULL
cna=as.data.frame(t(cna))
cna$SAMPLE_ID=row.names(cna)
cna=cna[adc_pts,c("SAMPLE_ID","ERBB2","BRCA2")]
cna=inner_join(meta,cna,by="SAMPLE_ID")
table(cna$ERBB2)
fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ ERBB2, data = cna)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

cox.test <- coxph(Surv(OS_MONTHS,OS_STATUS)~as.factor(ERBB2)+as.factor(STAGE_HIGHEST_RECORDED)+HR, data=cna) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)



pts=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_clinical_patient.txt")
sample=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_clinical_sample.txt")
meta=left_join(sample,pts,by="PATIENT_ID")%>%filter(CANCER_TYPE=="Breast Cancer",HER2=="Yes",SAMPLE_TYPE=="Metastasis")%>%as.data.frame()
meta$OS_STATUS[meta$OS_STATUS=="0:LIVING"]=0
meta$OS_STATUS[meta$OS_STATUS=="1:DECEASED"]=1
meta$OS_STATUS=as.numeric(meta$OS_STATUS)
agent=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_timeline_treatment.txt")
agent$duration=agent$STOP_DATE-agent$START_DATE
#msk_chord_2024,T-MD1
adc_pts=agent$PATIENT_ID[agent$AGENT%in%c("FAM-TRASTUZUMAB DERUXTECAN")&agent$START_DATE>0]
adc_pts=intersect(adc_pts,agent$PATIENT_ID[agent$AGENT!= c("PERTUZUMAB")])%>%unique()
adc_pts=intersect(meta$PATIENT_ID,adc_pts)
#tdm1_pts=intersect(agent$PATIENT_ID[agent$AGENT%in%c("TRASTUZUMAB")&agent$START_DATE<0&agent$duration>200],agent$PATIENT_ID[agent$AGENT%in%c("PERTUZUMAB")&agent$START_DATE<0&agent$duration>200])
meta=meta[meta$PATIENT_ID%in%adc_pts,]
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_cna.txt")%>%as.data.frame()
row.names(cna)=cna$Hugo_Symbol
cna$Hugo_Symbol=NULL
cna=as.data.frame(t(cna))
cna$SAMPLE_ID=row.names(cna)
cna=cna[adc_pts,c("SAMPLE_ID","ERBB2","BRCA2")]
cna=inner_join(meta,cna,by="SAMPLE_ID")
table(cna$ERBB2)
fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ ERBB2, data = cna)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

cox.test <- coxph(Surv(OS_MONTHS,OS_STATUS)~as.factor(ERBB2)+as.factor(STAGE_HIGHEST_RECORDED)+HR, data=cna) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)






