#Fig5b
library(forestplot)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
chisq.test(data$A01,data$lohhla)
variable=c("A01","A02","A03","A24","B07","B08","B27","B44")
data=data[!is.na(data$A01),]
results=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()
results=results[results$biomarker%in%variable,]
Total=results[,c("biomarker","whole_OR","Whole_LCI","Whole_UCI","whole_lr_p")]
Total$group="All"
Total$FDR=p.adjust(Total$whole_lr_p, method = "BH")
DHP=results[,c("biomarker","DHP_OR","DHP_LCI","DHP_UCI","DHP_lr_p")]
DHP$group="DHP"
DHP$FDR=p.adjust(DHP$DHP_lr_p, method = "BH")
TDM1=results[,c("biomarker","TDM1_OR","TDM1_LCI","TDM1_UCI","TDM1_lr_p")]
TDM1$group="T-DM1"
TDM1$FDR=p.adjust(TDM1$TDM1_lr_p, method = "BH")
colname=c("biomarker","OR","LCI","UCI","p","group","FDR")
colnames(Total)=colname;colnames(DHP)=colname;colnames(TDM1)=colname
df=rbind(DHP,TDM1)
df$biomarker=gsub("coding_mutation_", "",df$biomarker)
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
             clip = c(0,8),xlab = " OR with 95% CI") |> 
  fp_add_lines("black") |> 
  fp_add_header("HLA Supertype") |> 
  fp_set_style(box = c("#8491B4FF","#91D1C2FF") |> lapply(function(x) gpar(fill = x, col = "black")),
               default = gpar(vertices = TRUE)) |> 
  fp_set_zebra_style("#F5F9F9")


#Fig5c
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
#RNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
colnames(data)=gsub("Danaher-","",colnames(data))
colnames(data)=gsub("TIDE_","",colnames(data))
data1=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune.rds")
data1=data1[,c("patientID","Th2 cells","MHC.I_19272155")]
data=left_join(data1,data,by="patientID")%>%as.data.frame()
variable=c("B-cells","DC","Macrophages","T-cells","CD8-T-cells",     
           "Neutrophils","Cytotoxic-cells","Treg",
           "Mast-cells","NK-cells","CD45","TILs","Dysfunction",        
           "Exclusion","CAF","TAM_M2","Th2 cells","MHC.I_19272155",
           "TCR_clonality","BCR_clonality","FCGR3A","FCGR3B")
data[,variable]=scale(data[,variable])
data=data[!is.na(data$CAF),]
rna_results=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()
#DNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
variable=c("meanHED","lohhla","TCRA.tcell.fraction.adj","Neoantigen_DNA")
data$Neoantigen_DNA[is.na(data$Neoantigen_DNA)]=0
data[,c("meanHED","TCRA.tcell.fraction.adj")]=scale(data[,c("meanHED","TCRA.tcell.fraction.adj")])
data=data[!is.na(data$TCRA.tcell.fraction.adj),]
dna_results=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()

#image
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/image_metrics_PREDIX_HER2.rds")
varibale=c("Immune_Cell_prop","Distance_tumor_immune","Cell_Interaction")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin=clin[,c("patientID","Response","pCR","ER","Arm")]
data=left_join(data,clin,by="patientID")
image_results=Logistic_batch_adjER(data,"pCR","Arm",varibale,"ER")%>%as.data.frame()
# merge
results=rbind(rna_results,dna_results)%>%
  rbind(image_results)
write.csv(results, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure5/Figure5d.csv")

TDM1=results[,c("biomarker","TDM1_OR","TDM1_lr_p")]
TDM1$group="TDM1"
DHP=results[,c("biomarker","DHP_OR","DHP_lr_p")]
DHP$group="DHP"
colnames(TDM1)=c("Signature","OR","Pvalue","group")
colnames(DHP)=c("Signature","OR","Pvalue","group")
df=rbind(TDM1,DHP)
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
df$logOR=log(as.numeric(df$OR))
df$group=factor(df$group,levels = c("TDM1","DHP"))
unique(df$Signature)
df$Signature=factor(df$Signature,levels =rev(c("Neoantigen_DNA","meanHED","lohhla","MHC.I_19272155",
                                               "B-cells","DC","Macrophages","T-cells","CD8-T-cells","Cytotoxic-cells","NK-cells","CD45","TILs","Immune_Cell_prop",
                                               "TCR_clonality","BCR_clonality","TCRA.tcell.fraction.adj","FCGR3A","FCGR3B",
                                               "Dysfunction","Exclusion","CAF","Neutrophils","Mast-cells","Treg","TAM_M2","Th2 cells",
                                               "Distance_tumor_immune","Cell_Interaction")))
df=df[order(df$Signature,decreasing = F),]

baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=df,aes(Signature,logOR,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_modelbarplot(name="Treatment Arm")+
  labs(x="",y="lnOR")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))+
  coord_flip()

#5.5X6.5


#Fig5d
library(data.table)
library(tidyverse)
library(tableone)
library(forestploter)
require(openxlsx)
## continuous variable ##
# Distance_tumor_immune,Cell_Interaction,
# Th2 cells,FCGR3B Mast-cells Neutrophils  CAF
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
#RNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
colnames(data)=gsub("Danaher-","",colnames(data))
colnames(data)=gsub("TIDE_","",colnames(data))
data1=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune.rds")
data1=data1[,c("patientID","Th2 cells","MHC.I_19272155")]
data=left_join(data1,data,by="patientID")%>%as.data.frame()
variable=c("FCGR3B","CAF","Neutrophils","Mast-cells","Th2 cells")
data[,variable]=scale(data[,variable])
data=data[!is.na(data$CAF),]
immune=slice_metric(data,variable,20,80,5)%>%as.data.frame()
## re-run the batch logistic
selected_columns <- grep("_per_", colnames(immune), value = TRUE)
immune=immune[,c("pCR","Arm","ER",selected_columns)]
Interact_result=Logistic_batch_adjER(immune,"pCR","Arm",selected_columns,"ER")%>%as.data.frame()
immune$Arm=factor(immune$Arm,levels = c("DHP","T-DM1"))
res=Logistic_batch_continuous_subgroup(immune,selected_columns)%>%as.data.frame()
res_continuous=res[res$biomarker%in%c("FCGR3B_per_20","CAF_per_80","Mast-cells_per_35","Neutrophils_per_30","Th2 cells_per_60"),]
Interact_continuous=Interact_result[Interact_result$biomarker%in%c("FCGR3B_per_20","CAF_per_80","Mast-cells_per_35","Neutrophils_per_30","Th2 cells_per_60"),]
continuous=cbind(res_continuous,Interact_continuous)
#Digital no good results
#HLA supertype  Manually add
# merge
library(openxlsx);library(tableone)
list_of_datasets <- list("immune" = continuous)
write.xlsx(list_of_datasets, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Predictive_biomarker.xlsx")
table(genomic$Hypoxia_per_80,genomic$Arm)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
whole <- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$A01=="No",])
ShowRegTable(whole)
whole <- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$A01=="Yes",])
ShowRegTable(whole)
table(data$A01,data$Arm)
whole <- glm(as.numeric(pCR) ~ Arm+ER+A01+A01*Arm, family = "binomial", data = data)
ShowRegTable(whole)
# forest plot #
library(data.table)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure5/Subgroup_Forestplot.csv")
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

p <- forest(df[,c(1:3,7:9)],
            est = df$OR,
            lower = df$Low, 
            upper = df$High,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("DHP Better", "T-DM1 Better"),
            xlim = c(0, 4),
            ticks_at = c(0,0.5, 1, 2, 3,4,5),
            theme = tm)

# Print plot
plot(p)

#7X5




