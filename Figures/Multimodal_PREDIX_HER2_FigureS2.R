#############################################
##########Figure2a TMB ~ purity,ER#########
#############################################
library(ggpubr);library(ggplot2);library(data.table);library(tidyverse)
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
purity=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity_curated.csv"))
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic$patientID=as.integer(genomic$patientID)
data=left_join(genomic,purity,by="sampleID")%>%left_join(clin,by="patientID")
data=data[data$totalTMB*41.2>5,]
p1=ggscatter(data, x = "Purity", y = "TMB_uniform",
                add = "reg.line",  # Add regressin line
                xlim=c(0.1,0.75),ylim=c(0,20),
                add.params = list(color = "#0E4C92", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson", label.x = 0.3, label.y = 20) 


#4X5 p 
#############################################
######FigureS2b ERBB2/TP53_lollipop########
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
############FigureS2c-d forest plot############
#############################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
variable=c("coding_mutation_TP53","coding_mutation_PIK3CA","coding_mutation_GATA3","coding_mutation_ERBB2",
           "coding_mutation_TP53_oncokb","coding_mutation_PIK3CA_oncokb","coding_mutation_GATA3_oncokb")
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
df1=df[df$biomarker%in%c("TP53","PIK3CA","GATA3","ERBB2"),]
df1 |>
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

df2=df[df$biomarker%in%c("TP53_oncokb","PIK3CA_oncokb","GATA3_oncokb"),]
df2 |>
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

# 4x4

