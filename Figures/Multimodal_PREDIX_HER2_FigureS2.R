#############################################
##########Figure2A,B TMB ~ purity,ER#########
#############################################
library(ggpubr);library(ggplot2);library(data.table);library(tidyverse)
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
purity=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity_curated.csv"))
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic$patientID=as.integer(genomic$patientID)
data=left_join(genomic,purity,by="sampleID")%>%left_join(clin,by="patientID")
data=data[data$totalTMB*41.2>5,]
ggscatter(data, x = "Purity", y = "totalTMB",
                add = "reg.line",  # Add regressin line
                xlim=c(0.1,0.75),ylim=c(0,20),
                add.params = list(color = "#0E4C92", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson", label.x = 0.3, label.y = 20) 

#6X6 50%
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source(paste0(baseDir,"/Code/theme.R"))
d=data%>%select(c("ER","totalTMB"))
d <- reshape2::melt(d,id.vars=c("ER"))
figure_font_size=13
Fig <-
  ggplot(d,aes(x=ER,y=value,fill=ER))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_ER(name="ER status")+
  stat_compare_means(aes(group=ER),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.30)+
  labs(y="TMB(mut/MB)",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))+ 
  coord_cartesian(ylim = c(0,20))
Fig

#4X5 p 
#############################################
######FigureS2c ERBB2/TP53_lollipop########
#############################################
library(data.table);library(tidyverse)
maf=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_vaf0.05_curated.rds")
maf$vaf=maf$t_alt_count/maf$t_depth
vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                 "In_Frame_Ins", "Missense_Mutation")
maf=maf%>%filter(Variant_Classification%in%vc.nonSilent)
table(maf$ONCOGENIC)

a=maf[maf$Hugo_Symbol=="TP53"&maf$VARIANT_IN_ONCOKB==T,]
a=maf[maf$Hugo_Symbol=="ERBB2"&maf$VARIANT_IN_ONCOKB==T&maf$ONCOGENIC%in%c("Oncogenic","Likely Oncogenic"),]
table(maf$VARIANT_IN_ONCOKB,maf$ONCOGENIC)

maf=select(maf,c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification",'Reference_Allele',"Tumor_Seq_Allele2","Amino_acids"))
colnames(maf)=c("Sample_ID","Hugo_Symbol","Chromosome","Start_Position","End_Position","Mutation_Type",'Reference_Allele',"Variant_Allele","Protein_Change")
df=maf[maf$Hugo_Symbol%in%c("ERBB2"),]

write.table(df,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS2/ERBB2_TP53.txt",quote = F,row.names =F,sep="\t")

#HER2 mutation "UE-2971-1119-0" "UE-2971-1121-0" "UE-2971-1208-0" "UE-2971-1215-0" "UE-2971-1514-0" "UE-2971-1516-0" "UE-2971-1520-0"
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")

genomic$pCR2020[genomic$sampleID%in%c("UE-2971-1514-0","UE-2971-1516-0","UE-2971-1119-0",
                                      "UE-2971-1121-0","UE-2971-1208-0")]
table(genomic$coding_mutation_ERBB2)
#############################################
############FigureS2d forest plot############
#############################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
variable=c("coding_mutation_TP53","coding_mutation_PIK3CA","coding_mutation_ERBB2",
           "coding_mutation_TP53_oncokb","coding_mutation_PIK3CA_oncokb")
results=Logistic_batch_adjER(genomic,"pCR","Arm",variable,"ER")%>%as.data.frame()
colnames(results)
Total=results[,c("biomarker","whole_OR","Whole_LCI","Whole_UCI","whole_lr_p")]
Total$group="All"
DHP=results[,c("biomarker","DHP_OR","DHP_LCI","DHP_UCI","DHP_lr_p")]
DHP$group="DHP"
TDM1=results[,c("biomarker","TDM1_OR","TDM1_LCI","TDM1_UCI","TDM1_lr_p")]
TDM1$group="T-DM1"
colname=c("biomarker","OR","LCI","UCI","p","group")
colnames(Total)=colname;colnames(DHP)=colname;colnames(TDM1)=colname
df=rbind(Total,DHP,TDM1)
df$biomarker=gsub("coding_mutation_", "",df$biomarker)
#df$biomarker=gsub("_oncokb", "",df$biomarker)
biomarker=df$biomarker
df$OR=as.numeric(df$OR);df$LCI=as.numeric(df$LCI);df$UCI=as.numeric(df$UCI)
#
df |>
  group_by(group) |>
  forestplot(labeltext=biomarker,
             mean=OR,lower=LCI,upper=UCI,
             zero = 1,
             boxsize = .25, # We set the box size to better visualize the type
             line.margin = .1, # We need to add this to avoid crowding
             lty.ci = c(3),
             clip = c(0,4),xlab = " OR with 95% CI") |> 
  fp_add_lines("black") |> 
  fp_add_header("Mutation") |> 
  fp_set_style(box = c("#e5c06e", "#8491B4FF","#91D1C2FF") |> lapply(function(x) gpar(fill = x, col = "black")),
               default = gpar(vertices = TRUE)) |> 
  fp_set_zebra_style("#F5F9F9")


#TransNeo
meta=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=1)%>%as.data.frame()
meta=meta[meta$HER2.status=="POS"&meta$pCR.RD!="NA",]
maf=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=2)%>%as.data.frame()
maf$Tumor_Sample_Barcode=maf$Donor.ID;maf$Chromosome=maf$Chr;maf$Start_Position=maf$Start;maf$End_Position=maf$Start
maf$Reference_Allele=maf$Ref_Allele;maf$Tumor_Seq_Allele2=maf$Tumor_Alt_Allele
maf=maf[maf$Donor.ID%in%meta$Donor.ID,] 
write.table(maf,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/TransNeo_HER2_maf.txt",
            quote = F,sep = "\t",row.names = F)

maf=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/TransNeo_HER2_maf.oncokb_annotated.txt")
maf=maf[maf$VARIANT_IN_ONCOKB==T&maf$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")&maf$Hugo_Symbol=="TP53",]
meta=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=1)%>%as.data.frame()
meta=meta[meta$HER2.status=="POS"&meta$pCR.RD!="NA",]
meta$TP53="No"
meta$TP53[meta$Donor.ID%in%maf$Donor.ID]="Yes"
table(meta$TP53,meta$pCR.RD)
library(tableone)
meta$pCR=1
meta$pCR[meta$pCR.RD=="RD"]=0
res<- glm(pCR ~ TP53+ER.status, family = "binomial", data =meta)
ShowRegTable(res)
library(vcd);library("ggsci")
meta$pCR.RD=factor(meta$pCR.RD,levels = c("RD","pCR"))
meta$TP53=factor(meta$TP53,levels = c("Yes","No"))
mosaic( ~ ER.status + pCR.RD + TP53, data = meta,
        highlighting = "TP53", highlighting_fill = c("#BC3C29FF","#6F99ADFF"),
        direction = c("h","v","h"))


#############################################
############FigureS2e/f   KM plot############
#############################################
#DHP KM 
library("survival");library(survminer);library(tableone)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
table(genomic$coding_mutation_TP53_oncokb)
dhp=genomic[genomic$Arm=="DHP",]
dhp$TP53="Wild type"
dhp$TP53[dhp$coding_mutation_TP53_oncokb==1]="Mut"
dhp$TP53=factor(dhp$TP53,levels = c("Wild type","Mut"))
fit <- survfit(Surv(EFS.time,EFS.status) ~ TP53, data = dhp)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

table(dhp$coding_mutation_TP53_oncokb)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(TP53)+as.factor(Response)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=dhp) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)

#T-DM1 KM 
tdm1=genomic[genomic$Arm=="T-DM1",]
tdm1$TP53="Wild type"
tdm1$TP53[tdm1$coding_mutation_TP53_oncokb==1]="Mut"
tdm1$TP53=factor(tdm1$TP53,levels = c("Wild type","Mut"))
fit <- survfit(Surv(EFS.time,EFS.status) ~ TP53, data = tdm1)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))


table(tdm1$coding_mutation_TP53_oncokb)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(TP53)+as.factor(Response)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=tdm1) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)

# 6X4.5
########################################
#######TP53 mutation validation#########
# https://doi.org/10.1002/cam4.4652
########################################
#MSKCC -DHP
library("survival");library(survminer);library(tableone);library(data.table);library(tidyverse)
tumor=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_clinical_sample.txt")
#tumor=tumor[tumor$SAMPLE_TYPE=="Primary",]
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_clinical_patient.txt")
meta=left_join(meta,tumor,by="PATIENT_ID")
mut=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_mutations.oncokb_annotated.txt")
mut$patientID=substr(mut$Tumor_Sample_Barcode,1,9)
meta=meta[meta$PATIENT_ID%in%mut$patientID,]
mut=mut[mut$Hugo_Symbol=="TP53",]
mut=mut[mut$VARIANT_IN_ONCOKB==TRUE,]
mut=mut[mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")]
mut$vaf=mut$t_alt_count/(mut$t_alt_count+mut$t_ref_count)
range(mut$vaf)
mut=mut[mut$vaf>0.05,]
length(unique(mut$patientID))
meta$TP53="Wild type"
meta$TP53[meta$PATIENT_ID%in%mut$patientID]="Mut"
table(meta$TP53)

meta$TP53=factor(meta$TP53,levels = c("Wild type","Mut"))
meta$PFS.STATUS=substr(meta$PFS_STATUS,1,1)
meta$PFS.STATUS=as.numeric(meta$PFS.STATUS)
meta=meta%>%filter(!is.na(meta$PFS.STATUS)) 
#meta=meta[meta$DX_PSTAGE%in%c("IV"),]

fit <- survfit(Surv(PFS_MONTHS,PFS.STATUS) ~ TP53, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

table(meta$TP53)
cox.test <- coxph(Surv(PFS_MONTHS,PFS.STATUS)~as.factor(TP53)+as.factor(DX_PSTAGE), data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)
adc_pts=agent$PATIENT_ID[agent$AGENT%in%c("ADO-TRASTUZUMAB EMTANSINE")&agent$START_DATE>0]
adc_pts=intersect(adc_pts,agent$PATIENT_ID[agent$AGENT!= c("PERTUZUMAB")])%>%unique()
adc_pts=intersect(meta$PATIENT_ID,adc_pts)


#MSKCC T-DM1
pts=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_clinical_patient.txt")
sample=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_clinical_sample.txt")
meta=left_join(sample,pts,by="PATIENT_ID")%>%filter(CANCER_TYPE=="Breast Cancer",HER2=="Yes",SAMPLE_TYPE=="Primary",SOMATIC_STATUS=="Matched")%>%as.data.frame()
meta$OS_STATUS[meta$OS_STATUS=="0:LIVING"]=0
meta$OS_STATUS[meta$OS_STATUS=="1:DECEASED"]=1
meta$OS_STATUS=as.numeric(meta$OS_STATUS)
agent=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_timeline_treatment.txt")
agent$duration=agent$STOP_DATE-agent$START_DATE
adc_pts=agent$PATIENT_ID[agent$AGENT%in%c("ADO-TRASTUZUMAB EMTANSINE","FAM-TRASTUZUMAB DERUXTECAN")&agent$START_DATE>0]
adc_pts=intersect(adc_pts,agent$PATIENT_ID[agent$AGENT!= c("PERTUZUMAB")])%>%unique()
adc_pts=intersect(meta$PATIENT_ID,adc_pts)
meta=meta[meta$PATIENT_ID%in%adc_pts,]
mut=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/data_mutations.oncokb_annotated.txt")
sampleID=intersect(meta$SAMPLE_ID,mut$Tumor_Sample_Barcode)
meta=meta[meta$SAMPLE_ID%in%sampleID,];mut=mut[mut$Tumor_Sample_Barcode%in%sampleID,]
mut=mut[mut$VARIANT_IN_ONCOKB==TRUE&mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")&mut$Hugo_Symbol=="TP53",]

meta$TP53="Wild type"
meta$TP53[meta$SAMPLE_ID%in%mut$Tumor_Sample_Barcode]="Mut" 
table(meta$TP53)

meta$TP53=factor(meta$TP53,levels = c("Wild type","Mut"))
meta$OS.STATUS=substr(meta$OS_STATUS,1,1)%>%as.numeric()

fit <- survfit(Surv(OS_MONTHS,OS.STATUS) ~ TP53, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 100))
table(meta$TP53)
cox.test <- coxph(Surv(OS_MONTHS,OS.STATUS)~as.factor(TP53)+as.factor(HR)+CURRENT_AGE_DEID+CLINICAL_GROUP, data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)


tumor=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/msk_chord_2024/data_clinical_sample.txt")
tumor=tumor[tumor$SAMPLE_TYPE=="Primary"&tumor$PRIMARY_SITE=="Breast",]
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/msk_chord_2024/data_clinical_patient.txt")
meta=meta[meta$PATIENT_ID%in%tumor$PATIENT_ID,]
meta=meta[meta$GENDER=="Female"&meta$HER2=="Yes"&meta$STAGE_HIGHEST_RECORDED%in%c("Stage 1-3","Stage 4"),]
tumor=tumor[tumor$PATIENT_ID%in%meta$PATIENT_ID,]
meta=left_join(meta,tumor,by="PATIENT_ID")
trt=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/msk_chord_2024/data_timeline_treatment.txt")
trt=trt[trt$PATIENT_ID%in%meta$PATIENT_ID,]
trt$trt_duration=trt$STOP_DATE-trt$START_DATE
trt=trt[trt$trt_duration>100,]
dhp=intersect(trt$PATIENT_ID[trt$AGENT=="PERTUZUMAB"],trt$PATIENT_ID[trt$AGENT=="TRASTUZUMAB"])%>%unique()
dhp=c(dhp,unique(trt$PATIENT_ID[trt$AGENT=="HYALURONIDASE/PERTUZUMAB/TRASTUZUMAB"]))

meta=meta[meta$PATIENT_ID%in%dhp,]
tumor=tumor[tumor$PATIENT_ID%in%dhp,]
mut=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/msk_chord_2024/data_mutations.oncokb_annotated.txt")
mut=mut[mut$VARIANT_IN_ONCOKB==TRUE&mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic"),]
#mut=mut[mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")]
mut=mut[mut$Tumor_Sample_Barcode%in%tumor$SAMPLE_ID&mut$Hugo_Symbol=="TP53",]
tumor=tumor[tumor$SAMPLE_ID%in%mut$Tumor_Sample_Barcode,]

meta$TP53="Wild type"
meta$TP53[meta$PATIENT_ID%in%tumor$PATIENT_ID]="Mut" 
table(meta$TP53)

meta$TP53=factor(meta$TP53,levels = c("Wild type","Mut"))
meta$OS.STATUS=substr(meta$OS_STATUS,1,1)%>%as.numeric()

fit <- survfit(Surv(OS_MONTHS,OS.STATUS) ~ TP53, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(meta$TP53)
cox.test <- coxph(Surv(OS_MONTHS,OS.STATUS)~as.factor(TP53)+as.factor(HR)+CURRENT_AGE_DEID+CLINICAL_GROUP, data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)
