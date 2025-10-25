#Fig.2a
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(data.table)
#Mutation
source("E:/Projects/PREDIX_HER2/Multimodal/Code/maf2oncoprint.R")
mutect=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_curated.txt")
maf=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_vaf0.05_curated.rds")
drivers=read_excel("E:/Projects/PREDIX_HER2/CUTseq/data/genelist/NIHMS68344-supplement-Supplementary_Tables_1-21/nature17676-s3/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx",sheet=4)
drivers_mut=drivers$Gene[drivers$Mutation_Type%in%c("Del","Ins","Rearrangement","Sub")]
drivers_cna=drivers$Gene[drivers$Mutation_Type%in%c("CopyNumber")]%>%unique()
table(drivers$Mutation_Type)
str(maf)
# mannual driver
HER_pathway=c("EGFR","ERBB2","ERBB3","ERBB4") # mut
TP53_pathway=c("ATM","TP53") # mut
PIK3_AKT_pathway=c("PIK3CA","PIK3R1","PTEN","AKT1") # PTEN del
MAPK_ERK_pathway=c("MAPK","MAP2K","MAP3K1","MAP2K4","NF1") # NF1 del
CDK_RB_pathway=c("CCND1","CDK4","CDK6","RB1") # CCND1 and CDK4/6 AMP  RB1 Del
###
#table(maf$Variant_Classification)
counts <- table(maf$Hugo_Symbol)
gene <- names(counts[counts >5])
drivers_mut=intersect(gene,unique(drivers_mut))
drivers_mut=c(drivers_mut,HER_pathway,TP53_pathway,PIK3_AKT_pathway,MAPK_ERK_pathway,CDK_RB_pathway)%>%unique()
drivers_mut=intersect(drivers_mut,maf$Hugo_Symbol)
#drivers=intersect(drivers,unique(rownames(data)))
##SCNA
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_thresholded.by_genes.txt")%>%filter(`Gene Symbol`%in% drivers_mut)%>%as.data.frame()
rownames(CNA)=CNA$`Gene Symbol`
CNA=as.data.frame(t(CNA))
CNA=CNA[4:nrow(CNA),]
sampleID=row.names(CNA)%>%as.character()
CNA=sapply(CNA, as.numeric)%>%as.data.frame()
CNA[CNA=="-2"]="Deletion"
CNA[CNA=="-1"]="Neutral"
CNA[CNA=="0"]="Neutral"
CNA[CNA=="1"]="Neutral"
CNA[CNA=="2"]="Amplification"
CNA$sampleID=sampleID
pid=data.frame(sampleID=unique(mutect$Tumor_Sample_Barcode))
CNA=inner_join(pid,CNA,by="sampleID")
CNA$patientID=NULL
library(reshape)
CNA_long <- melt(CNA,id=c("sampleID"))%>%as.data.frame()
CNA_long=CNA_long[CNA_long$value!="Neutral",]
colnames(CNA_long)=c("Tumor_Sample_Barcode",'Hugo_Symbol',"Variant_Classification")
CNA_long=CNA_long[CNA_long$Hugo_Symbol%in%drivers_mut,]
## Merge mutation and CNA
maf=rbind(maf,CNA_long,fill=TRUE)
maf[is.na(maf)]=1
tail(maf)
maf=as.data.frame(maf)
## maf2concoprint
setDT(maf)
data=maf2oncoprint(maf)
## TMB and genomic variables ##
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=genomic[genomic$totalTMB!=0,]
col <- c("Missense" = "#00A087FF","Truncating" = "#ac9141","In.frame"="#F39B7FFF","Silent" = "#916ba6",
         "Amplification"="#BC102B","Deletion"="#5385AC") #"Gain"="#DF989E","Loss"="#AEC4D6",

alter_fun <- list(
  background = function(x,y,w,h){
    grid::grid.rect(x,y,w - grid::unit(0.5,"mm"), h - grid::unit(0.5,"mm"),gp=grid::gpar(fill = "#CCCCCC",col=NA))
  },
  Missense = function(x,y,w,h){
    grid::grid.rect(x,y,w - grid::unit(0.5,"mm"),h - grid::unit(0.5, "mm"),gp=grid::gpar(fill=col["Missense"],col=NA))
  },
  Truncating = function(x,y,w,h){
    grid::grid.rect(x,y,w - grid::unit(0.5,"mm"),h - grid::unit(0.5,"mm"),gp = grid::gpar(fill=col["Truncating"],col=NA))
  },
  In.frame = function(x,y,w,h){
    grid::grid.rect(x,y,w - grid::unit(0.5,"mm"),h - grid::unit(0.5,"mm"),gp = grid::gpar(fill=col["In.frame"],col=NA))
  },
  Silent = function(x,y,w,h){
    grid::grid.rect(x,y,w - grid::unit(0.5,"mm"),h - grid::unit(0.5,"mm"),gp = grid::gpar(fill=col["Silent"],col=NA))
  },
  Amplification = function(x, y, w, h) {grid.rect(x, y, w-unit(2, "pt"), h*0.33,gp = gpar(fill = col["Amplification"], col = NA))
  },
  Deletion = function(x, y, w, h) {grid.rect(x, y, w-unit(2, "pt"), h*0.33,gp = gpar(fill = col["Deletion"], col = NA))
  }
)
clinical_anno=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clinical_anno=inner_join(clinical_anno,genomic,by='patientID')
clinical_anno$log2TMB=as.numeric(log2(clinical_anno$totalTMB))
clinical_anno=clinical_anno[order(clinical_anno$TREAT,-clinical_anno$log2TMB),]
clinical_anno=clinical_anno[clinical_anno$patientID%in%unique(substr(colnames(data),9,12)),]
sample_order=clinical_anno$sampleID
##heatmap
gene_order=c("ERBB2","TP53","PIK3CA","MYC","GNAS","ZNF217","CCND1","NCOR1","MAP2K4","GATA3","FOXA1","CDH1","NF1",
             "ATM","MAP3K1","KRAS","EGFR","BCOR","ARID1A","XBP1","CCNE1","RUNX1","SPEN","CUX1","TBX3","BRCA2","CIC",
             "SMARCA4","ASXL1","MDM2","ATR","CBLB","KDM6A","ZFP36L1","CNOT3","CREBBP","SETD2","TET2","APC","STAG2",
             "USP9X","PTEN","BRAF","AKT1")
HER_pathway=c("EGFR","ERBB2","ERBB3","ERBB4") # mut
TP53_pathway=c("ATM","TP53") # mut
PIK3_AKT_pathway=c("PIK3CA","PIK3R1","PTEN","AKT1") # PTEN del
MAPK_ERK_pathway=c("MAPK","MAP2K","MAP3K1","MAP2K4","NF1") # NF1 del
CDK_RB_pathway=c("CCND1","CDK4","CDK6","RB1") # CCND1 and CDK4/6 AMP  RB1 Del
drivers_mut=intersect(drivers_mut,row.names(data))
row_annot=data.frame(gene=drivers_mut)
row_annot$HER_pathway="No"
row_annot$HER_pathway[row_annot$gene%in%HER_pathway]="Yes"
row_annot$TP53_pathway="No"
row_annot$TP53_pathway[row_annot$gene%in%TP53_pathway]="Yes"
row_annot$PIK3_AKT_pathway="No"
row_annot$PIK3_AKT_pathway[row_annot$gene%in%PIK3_AKT_pathway]="Yes"
row_annot$MAPK_ERK_pathway="No"
row_annot$MAPK_ERK_pathway[row_annot$gene%in%MAPK_ERK_pathway]="Yes"
row_annot$CDK_RB_pathway="No"
row_annot$CDK_RB_pathway[row_annot$gene%in%CDK_RB_pathway]="Yes"
##left annotation
##plot heatmap
ht=oncoPrint(data[drivers_mut,sample_order],alter_fun = alter_fun, col = col,column_order = sample_order,alter_fun_is_vectorized = FALSE,
             pct_side = "right", row_names_side = "left",row_names_gp = grid::gpar(fontsize = 12),
             left_annotation = rowAnnotation(TP53_pathway =as.factor(row_annot$TP53_pathway),
                                             HER_pathway =as.factor(row_annot$HER_pathway), 
                                             PIK3_AKT_pathway=as.factor(row_annot$PIK3_AKT_pathway),
                                             MAPK_ERK_pathway=as.factor(row_annot$MAPK_ERK_pathway),
                                             CDK_RB_pathway=as.factor(row_annot$CDK_RB_pathway),
                                             col=list(TP53_pathway=c("Yes"="#525252","No"="#f0f0f0"),
                                                      HER_pathway=c("Yes"="#525252","No"="#f0f0f0"),PIK3_AKT_pathway=c("Yes"="#525252","No"="#f0f0f0"),
                                                      MAPK_ERK_pathway=c("Yes"="#525252","No"="#f0f0f0"),CDK_RB_pathway=c("Yes"="#525252","No"="#f0f0f0")),
                                             border=TRUE,simple_anno_size = unit(2, "mm")),
             top_annotation = HeatmapAnnotation(simple_anno_size = unit(2, "mm"),
                                                Uniform_TMB= anno_barplot(clinical_anno$log2TMB,border = F,bar_width = 0.4),
                                                Treatment=clinical_anno$Arm,
                                                Response=clinical_anno$Response,
                                                Recurrence=clinical_anno$EFS.status,
                                                HER2_IHC=clinical_anno$HER2neu1,
                                                HR=clinical_anno$ERPRdic,
                                                col=list(Treatment=c("T-DM1"="#91D1C2FF","DHP"="#8491B4FF"),
                                                         Response=c("pCR"="#f0f0f0","RD"="#525252"),
                                                         Recurrence=c("0"="#f0f0f0","1"="#525252"),
                                                         HR=c("ERandPR-"="#f0f0f0","ERorPR+"="#525252"),
                                                         HER2_IHC=c("2+"="#f0f0f0","3+"="#525252"),
                                                         na_col="grey",border=TRUE)))
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "bottom")
# 13X 7 70%
#%v%Heatmap(CIN[,substr(sample_order,9,12)], col = colorRamp2(c(-1,-0.5,0,0.5,3),c("#5385AC", "#AEC4D6", "#EAEAEA","#DF989E","#BC102B")),show_column_names = FALSE,
#          name = "CINSignature", height = unit(1.5, "cm"),row_names_gp = grid::gpar(fontsize = 9))%v%
#  Heatmap(COSMIC[,substr(sample_order,9,12)], col = colorRamp2(c(-1,-0.5,0,0.5,3),c("#5385AC", "#AEC4D6", "#EAEAEA","#DF989E","#BC102B")), 
#          name = "COSMIC.signature",show_column_names = FALSE, height = unit(2, "cm"),row_names_gp = grid::gpar(fontsize = 7.5))
# export source data
df=data[drivers_mut,sample_order]%>%t()%>%as.data.frame()
df$sampleID=row.names(df)
clin=clinical_anno[,c("sampleID","log2TMB","Arm","Response","RFS.status","HER2neu1","ER")]
colnames(clin)[1]="sampleID"
df=left_join(clin,df,by="sampleID")
write.table(df,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure2/source_data_figure2a.txt",quote = F,row.names =F,sep="\t")

#Fig.2b
library(vcd);library("ggsci");library(tidyverse);library(data.table)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
genomic=left_join(genomic,clin,by="patientID") # genomic=genomic[genomic$switch=="No",]
genomic$Response=factor(genomic$Response,levels = c("RD","pCR"))
genomic$TP53[genomic$coding_mutation_TP53_oncokb==1]="Mutation"
genomic$TP53[genomic$coding_mutation_TP53_oncokb==0]="Wild type"


ftable(genomic$Arm,genomic$ER,genomic$TP53,genomic$Response)
ftable(genomic$Arm)
mosaic( ~ Arm + ER + Response + TP53, data = genomic,
        highlighting = "TP53", highlighting_fill = c("#BC3C29FF","#6F99ADFF"),
        direction = c("v","h","v","h"))

# 5X5 55%
DHP=genomic[genomic$Arm=="DHP",]
TDM1=genomic[genomic$Arm=="T-DM1",]

table(DHP$coding_mutation_TP53_oncokb,DHP$Response,DHP$ER)
table(TDM1$coding_mutation_TP53_oncokb,TDM1$Response,TDM1$ER)

interaction_2<- glm(pCR ~ coding_mutation_TP53_oncokb*Arm+coding_mutation_TP53_oncokb+Arm+ER, family = "binomial", data =genomic)
library(tableone)
ShowRegTable(interaction_2)


interaction_2<- glm(pCR ~ coding_mutation_TP53_oncokb+ER, family = "binomial", data =genomic[genomic$Arm=="T-DM1",])
library(tableone)
ShowRegTable(interaction_2)
interaction_2<- glm(pCR ~ coding_mutation_TP53_oncokb+ER, family = "binomial", data =genomic[genomic$Arm=="DHP",])
library(tableone)
ShowRegTable(interaction_2)
#Fig.2c
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
DHP=results[,c("biomarker","group","DHP_OR","DHP_LCI","DHP_UCI","DHP_lr_p")]
DHP$DHP_OR=as.numeric(DHP$DHP_OR)
DHP$DHP_lr_p=as.numeric(DHP$DHP_lr_p)
DHP$FDR=p.adjust(DHP$DHP_lr_p, method = "BH")
DHP$lnOR=DHP$DHP_OR%>%log(base=exp(1))
DHP$log10FDR=-log10(DHP$FDR)

library(ggplot2)
library(ggrepel)
# FDR and color
thr_5FDR <- -log10(0.05)
color_map <- c("amp" = "brown", "del" = "steelblue")
# 绘图
p1 <- ggplot(DHP, aes(x = lnOR, y = log10FDR, fill= group)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.3,
             position = position_jitter(width = 0.05, height = 0.02)) +  # 外框点
  scale_fill_manual(values = color_map) +
  geom_hline(yintercept = thr_5FDR, color = "brown", linetype = "solid") +
  geom_vline(xintercept = 0, color = "gray40", linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(DHP, FDR < 0.25),
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
#--------------
#Fig.2d
#--------------
TDM1=results[,c("biomarker","group","TDM1_OR","TDM1_LCI","TDM1_UCI","TDM1_lr_p")]
TDM1$TDM1_OR=as.numeric(TDM1$TDM1_OR)
TDM1$TDM1_lr_p=as.numeric(TDM1$TDM1_lr_p)
TDM1$FDR=p.adjust(TDM1$TDM1_lr_p, method = "BH")
TDM1$lnOR=TDM1$TDM1_OR%>%log(base=exp(1))
TDM1$log10FDR=-log10(TDM1$FDR)

library(ggplot2)
library(ggrepel)
# FDR and color
thr_5FDR <- -log10(0.05)
color_map <- c("amp" = "brown", "del" = "steelblue")
# 绘图
p2 <- ggplot(TDM1, aes(x = lnOR, y = log10FDR, fill= group)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.3,
             position = position_jitter(width = 0.05, height = 0.02)) +  # 外框点
  scale_fill_manual(values = color_map) +
  geom_hline(yintercept = thr_5FDR, color = "brown", linetype = "solid") +
  geom_vline(xintercept = 0, color = "gray40", linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(TDM1, FDR < 0.25),
                           aes(label = biomarker),
                           size = 4, color = "black") +
  annotate("text", x = 2, y = thr_5FDR + 0.05, label = "5% FDR", color = "brown", hjust = 1) +
  annotate("text", x = -3, y = 0.2, label = "RD associated", color = "black", hjust = 0, size = 4) +
  annotate("text", x = 3, y = 0.2, label = "pCR associated", color = "black", hjust = 1, size = 4) +
  
  labs(x = "ln(ORadjusted)", y = expression(-log[10](FDR))) +
  xlim(-3, 3) +
  ylim(0, 2) +
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank())

p2


library(ggpubr)
ggarrange(p1,p2)
#8X3.5

#--------------
#Fig.2e
#--------------
##################################
########P value genomic metrics###
##################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
data=left_join(genomic,clin,by="patientID")
variable0=c("COSMIC.Signature.2","COSMIC.Signature.3","COSMIC.Signature.6",
            "COSMIC.Signature.7","COSMIC.Signature.10","COSMIC.Signature.13")
data[,variable0]=data[,variable0] %>% mutate(across(where(is.numeric), scale))
df=data[!is.na(data$COSMIC.Signature.2),]
results0=Logistic_batch_adjER(df,"pCR","Arm",variable0,"ER")%>%as.data.frame()

variable=c("TMB_uniform","TMB_clone","CNV_burden","LOH_Del_burden","HRD")
data[,variable]=data[,variable] %>% mutate(across(where(is.numeric), scale))
df=data
results=Logistic_batch_adjER(df,"pCR","Arm",variable,"ER")%>%as.data.frame()
results=rbind(results,results0)

fig2e=results
TDM1=results[,c("biomarker","TDM1_OR","TDM1_lr_p")]
TDM1$group="TDM1"
DHP=results[,c("biomarker","DHP_OR","DHP_lr_p")]
DHP$group="DHP"
colnames(TDM1)=c("Genomic","OR","Pvalue","group")
colnames(DHP)=c("Genomic","OR","Pvalue","group")
df=rbind(TDM1,DHP)
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
df$logOR=log(as.numeric(df$OR))
df$group=factor(df$group,levels = c("TDM1","DHP"))
unique(df$Genomic)
df$Genomic=factor(df$Genomic,levels =c("HRD","LOH_Del_burden","CNV_burden","COSMIC.Signature.13","COSMIC.Signature.10","COSMIC.Signature.7","COSMIC.Signature.6",
                                       "COSMIC.Signature.3","COSMIC.Signature.2","TMB_clone","TMB_uniform"))
df=df[order(df$Genomic),]

baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=df,aes(Genomic,logOR,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_modelbarplot(name="Treatment Arm")+
  labs(x="Genomic Metrics",y="lnOR")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,1), "lines"))+
  coord_flip()

#7X4


#Fig.2f
##################################
#########Subgroup forest##########
##################################
library(data.table)
library(tidyverse)
library(tableone)
library(forestploter)
require(openxlsx)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
## binary variable ##
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")%>%as.data.frame()
CNA=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/CUTseq_gene_baseline_complemental.rds")
mutation=genomic[,c("patientID","coding_mutation_TP53_oncokb")]
CNA=CNA[,c("BRCA2","FGFR1","NCOR1","ERBB2","PIK3CA","MAP3K1")]
colnames(CNA)=paste0("CNA_",colnames(CNA))
CNA$patientID=row.names(CNA)%>%as.integer()
data=left_join(clin,mutation,by="patientID")%>%left_join(CNA,by="patientID")%>%as.data.frame()
data$CNA_BRCA2_Del[data$CNA_BRCA2%in%c(0,1,2)]="0"
data$CNA_BRCA2_Del[data$CNA_BRCA2%in%c(-2,-1)]="1"
data$CNA_ERBB2_Amp[data$CNA_ERBB2%in%c(0,-1,-2)]="0"
data$CNA_ERBB2_Amp[data$CNA_ERBB2%in%c(2,1)]="1"

bin_var=c("coding_mutation_TP53_oncokb","CNA_BRCA2_Del","CNA_ERBB2_Amp")
data[,bin_var]=lapply(as.data.frame(data[,bin_var]),function(x) as.factor(x))
data$Arm=factor(data$Arm,levels = c("DHP","T-DM1"))
res=Logistic_batch_bin_subgroup(data,bin_var)%>%as.data.frame()

df=data[!is.na(data$coding_mutation_TP53_oncokb),]
variable=c("coding_mutation_TP53_oncokb")
mut=Logistic_batch_adjER(df,"pCR","Arm",variable,"ER")%>%as.data.frame()

df=data[!is.na(data$CNA_BRCA2_Del),]
variable=c("CNA_BRCA2_Del","CNA_ERBB2_Amp")
cna=Logistic_batch_adjER(df,"pCR","Arm",variable,"ER")%>%as.data.frame()
res_bin=rbind(mut,cna)

all.equal(res$biomarker,res_bin$biomarker)
bin=cbind(res,res_bin)
## continuous variable ##
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")%>%as.data.frame()
genomic=inner_join(genomic,clin,by="patientID")%>%as.data.frame()
variable=c("LOH_Del_burden","COSMIC.Signature.13")
genomic=slice_metric(genomic,variable,20,80,5)
## re-run the batch logistic
selected_columns <- grep("_per_", colnames(genomic), value = TRUE)
Interact_result=Logistic_batch_adjER(genomic,"pCR","Arm",selected_columns ,"ER")%>%as.data.frame()
genomic$Arm=factor(genomic$Arm,levels = c("DHP","T-DM1"))
res=Logistic_batch_continuous_subgroup(genomic,selected_columns)%>%as.data.frame()
res_continuous=res[res$biomarker%in%c("LOH_Del_burden_per_60","COSMIC.Signature.13_per_55"),]
Interact_continuous=Interact_result[Interact_result$biomarker%in%c("LOH_Del_burden_per_60","COSMIC.Signature.13_per_55"),]
continuous=cbind(res_continuous,Interact_continuous)
# merge
library(openxlsx)
list_of_datasets <- list("bin" = bin, "continuous" = continuous)
write.xlsx(list_of_datasets, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure2/Predictive_biomarker.xlsx")
# forest plot #
library(data.table)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure1/Subgroup_Forestplot.csv")
df$DHP <- ifelse(is.na(df$DHP), "", df$DHP)
df$`T-DM1` <- ifelse(is.na(df$`T-DM1`), "", df$`T-DM1`)
df$`P for interaction` <- ifelse(is.na(df$`P for interaction`), "", df$`P for interaction`)
df$` ` <- paste(rep(" ", 20), collapse = " ")
a=df$`OR (95% CI)`;b=df$`P for interaction`
df$`OR (95% CI)`=NULL;df$`P for interaction`=NULL
df$`OR (95% CI)`=a
df$`P for interaction`=b

tm <- forest_theme(base_size = 10,
                   refline_col = "black",
                   arrow_type = "closed",
                   footnote_col = "blue")

p <- forest(df[,c(1:3,8:10)],
            est = df$OR,
            lower = df$Low, 
            upper = df$High,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("DHP Better", "T-DM1 Better"),
            xlim = c(0, 4),
            ticks_at = c(0,0.5, 1, 2, 3),
            theme = tm)

# Print plot
plot(p)

#7X7


# export data

library(openxlsx)
data_list=list("fig2b"=fig2b%>%as.data.frame(),
               "fig2c"=fig2c%>%as.data.frame(),
               "fig2d"=fig2d%>%as.data.frame())
openxlsx::write.xlsx(data_list,file='E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Figure2.xlsx')



