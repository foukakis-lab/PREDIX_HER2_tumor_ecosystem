# fig7a
#########feature correlation##########
library(data.table);library(ggplot2);library(ggpubr);library(tidyverse);library(ComplexHeatmap)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
#data=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
data=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2_withNA.txt")%>%as.data.frame()
data$WSI_Immune_Cell_prop=as.numeric(data$WSI_Immune_Cell_prop)
data$WSI_Distance_tumor_immune=as.numeric(data$WSI_Distance_tumor_immune)
data$WSI_Cell_Interaction=as.numeric(data$WSI_Cell_Interaction)
data$RNA_sspbc.subtype=NULL
rna_cols <- grep("^RNA_", colnames(data), value = TRUE)
dna_cols <- grep("^DNA_", colnames(data), value = TRUE)
prot_cols <- grep("^Prot_", colnames(data), value = TRUE)
wsi_cols <- grep("^WSI_", colnames(data), value = TRUE)
data$Arm=data$Clin_Arm

# univariate association---RNA
df=data[,c(rna_cols,"pCR","Arm")]
df <- df[complete.cases(df), ]
results_rna=Logistic_batch_uni(df,"pCR","Arm",rna_cols)%>%as.data.frame()
results=results_rna
sig_vars <- data.frame(
  biomarker = results$biomarker,
  whole_sig = results$whole_lr_p < 0.05,
  DHP_sig   = results$DHP_lr_p < 0.05,
  TDM1_sig  = results$TDM1_lr_p < 0.05
) %>%
  filter(whole_sig | DHP_sig | TDM1_sig) %>%
  pull(biomarker)
results_rna=results[results$biomarker%in%sig_vars,]
# univariate association---DNA
df=data[,c(dna_cols,"pCR","Arm")]
df$DNA_Neoantigen_DNA[is.na(df$DNA_Neoantigen_DNA)]=0
df <- df[!is.na(df$DNA_ERBB2_CNA), ]
results_dna=Logistic_batch_uni(df,"pCR","Arm",dna_cols)%>%as.data.frame()
results=results_dna
sig_vars <- data.frame(
  biomarker = results$biomarker,
  whole_sig = results$whole_lr_p < 0.05,
  DHP_sig   = results$DHP_lr_p < 0.05,
  TDM1_sig  = results$TDM1_lr_p < 0.05
) %>%
  filter(whole_sig | DHP_sig | TDM1_sig) %>%
  pull(biomarker)
results_dna=results[results$biomarker%in%sig_vars,]
# univariate association---Proteomics
df=data[,c(prot_cols,"pCR","Arm")]
df <- df[complete.cases(df), ]
results_prot=Logistic_batch_uni(df,"pCR","Arm",prot_cols)%>%as.data.frame()
results=results_prot
sig_vars <- data.frame(
  biomarker = results$biomarker,
  whole_sig = results$whole_lr_p < 0.05,
  DHP_sig   = results$DHP_lr_p < 0.05,
  TDM1_sig  = results$TDM1_lr_p < 0.05
) %>%
  filter(whole_sig | DHP_sig | TDM1_sig) %>%
  pull(biomarker)
results_prot=results[results$biomarker%in%sig_vars,]
# univariate association---wsi
df=data[,c(wsi_cols,"pCR","Arm")]
df <- df[complete.cases(df), ]
results_wsi=Logistic_batch_uni(df,"pCR","Arm",wsi_cols)%>%as.data.frame()
results=results_wsi
sig_vars <- data.frame(
  biomarker = results$biomarker,
  whole_sig = results$whole_lr_p < 0.05,
  DHP_sig   = results$DHP_lr_p < 0.05,
  TDM1_sig  = results$TDM1_lr_p < 0.05
) %>%
  filter(whole_sig | DHP_sig | TDM1_sig) %>%
  pull(biomarker)
results_wsi=results[results$biomarker%in%sig_vars,]

results=rbind(results_rna,results_dna)%>%rbind(results_prot)%>%rbind(results_wsi)
#############
results <- results %>%
  mutate(All = case_when(
    whole_OR.OR > 1 & whole_lr_p < 0.05 ~ "Positive",
    whole_OR.OR < 1 & whole_lr_p < 0.05 ~ "Negative",
    TRUE ~ "NS"
  ))

results <- results %>%
  mutate(DHP = case_when(
    DHP_OR.OR > 1 & DHP_lr_p < 0.05 ~ "Positive",
    DHP_OR.OR < 1 & DHP_lr_p < 0.05 ~ "Negative",
    TRUE ~ "NS"
  ))

results <- results %>%
  mutate(TDM1 = case_when(
    TDM1_OR.OR > 1 & TDM1_lr_p < 0.05 ~ "Positive",
    TDM1_OR.OR < 1 & TDM1_lr_p < 0.05 ~ "Negative",
    TRUE ~ "NS"
  ))
results$metrics=substr(results$biomarker,1,3)
results$metrics=factor(results$metrics,levels = c("DNA","RNA","Pro","WSI"))
ha=HeatmapAnnotation(Feature=results$metrics,
                     All=results$All,
                     DHP=results$DHP,
                     TDM1=results$TDM1,
                     col=list(Feature=c("DNA"="#f2a104","RNA"="#72a2c0","Pro"="#9467bd","WSI"="#00743f"),
                              All=c("Positive"="#00A087FF","Negative"="#fdb462","NS"="#BEBEBE"),
                              DHP=c("Positive"="#00A087FF","Negative"="#fdb462","NS"="#BEBEBE"),
                              TDM1=c("Positive"="#00A087FF","Negative"="#fdb462","NS"="#BEBEBE")),
                     na_col="#808080",show_annotation_name = FALSE,
                     annotation_height = unit(c(0.5,0.5,0.5,0.5), "mm"),
                     border=F)

mat=data[,results$biomarker]
mat$DNA_coding_mutation_TP53_oncokb[mat$DNA_coding_mutation_TP53_oncokb==FALSE]=0
mat$DNA_coding_mutation_TP53_oncokb[mat$DNA_coding_mutation_TP53_oncokb==TRUE]=1
cor_matrix <- cor(mat, method = "pearson", use = "complete.obs")
min(cor_matrix);max(cor_matrix)

library(circlize)
blue2red <- colorRampPalette(c("midnightblue","#003BC0","#6689d9","white",
                               "#eb8694","#D80D29","firebrick4"))
breaks <- seq(-1, 1, length.out = 7)
colors <- blue2red(7)
col_fun <- colorRamp2(breaks, colors)

p=Heatmap(
  cor_matrix, 
  name = "Correlation", 
  col = col_fun , 
  top_annotation = ha,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = TRUE, 
  show_column_names = TRUE, 
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  show_heatmap_legend = TRUE
)
p
pdf("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Fig7a.pdf", width = 9, height = 8)
draw(p)  
dev.off()

# fig7b-d
##########All##########
library(data.table);library(ggplot2);library(Polychrome)
Int=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Both_Integrate/all_no_norm_selfeats_10_run4_ROC_values.csv")
Clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Both_Treat_clinical/all_no_norm_selfeats_6_clinical_ROC_values.csv")
DNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Both_Treats_DNA/all_no_norm_selfeats_10_DNA_ROC_values.csv")
RNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Both_Treats_RNA/all_no_norm_selfeats_10_RNA_ROC_values.csv")
Prot=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Both_Proteomics/all_no_norm_selfeats_10_Proteomics_ROC_values.csv")
WSI=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Both_Treat_Image/all_no_norm_selfeats_3_Image_ROC_values.csv")

rename <- c("fpr","tpr")

names(Int)[2:3] <- rename
names(Clin)[2:3] <- rename
names(DNA)[2:3] <- rename
names(RNA)[2:3] <- rename
names(Prot)[2:3] <- rename
names(WSI)[2:3] <- rename

Clin[, data := "Clinic Only - ROC AUC: 0.67 ± 0.11"]
DNA[, data := "Genomic Only - ROC AUC: 0.71 ± 0.11"]
RNA[, data := "Transcriptomic Only - ROC AUC: 0.83 ± 0.1"]
Prot[, data := "Proteomic Only - ROC AUC: 0.75 ± 0.1"]
WSI[, data := "Digital Pathology Only - ROC AUC: 0.57 ± 0.11"]
Int[, data := "Integrated - ROC AUC: 0.88 ± 0.08"]


All <- rbind(Clin,DNA,RNA,Prot,WSI,Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Clinic Only - ROC AUC: 0.67 ± 0.11","Genomic Only - ROC AUC: 0.71 ± 0.11",
                                    "Transcriptomic Only - ROC AUC: 0.83 ± 0.1","Proteomic Only - ROC AUC: 0.75 ± 0.1",
                                    "Digital Pathology Only - ROC AUC: 0.57 ± 0.11","Integrated - ROC AUC: 0.88 ± 0.08"))

ROC_all <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +     
  scale_colour_manual(values = c("#434279", "#f2a104","#72a2c0","#9467bd","#00743f","#d62728")) +# protein "#9467bd"
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_all

ggsave(ROC_all, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/All_roc_curves.pdf", width = 6.2, height = 6.2)


##########T-DM1##########
library(data.table);library(ggplot2);library(Polychrome)
Int=fread("E:Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_1_Integrated/all_no_norm_TREAT_1_selfeats_10_run4_ROC_values.csv")
Clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_1_clinical/all_no_norm_TREAT_1_selfeats_6_clinical_ROC_values.csv")
DNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_1_DNA/all_no_norm_TREAT_1_selfeats_10_DNA_ROC_values.csv")
RNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_1_RNA/all_no_norm_TREAT_1_selfeats_10_RNA_ROC_values.csv")
Prot=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_1_Proteomics/all_no_norm_TREAT_1_selfeats_10_Proteomics_ROC_values.csv")
WSI=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_1_Image/all_no_norm_TREAT_1_selfeats_3_Image_ROC_values.csv")

rename <- c("fpr","tpr")

names(Int)[2:3] <- rename
names(Clin)[2:3] <- rename
names(DNA)[2:3] <- rename
names(RNA)[2:3] <- rename
names(Prot)[2:3] <- rename
names(WSI)[2:3] <- rename

Clin[, data := "Clinic Only - ROC AUC: 0.61 ± 0.16"]
DNA[, data := "Genomic Only - ROC AUC: 0.79 ± 0.15"]
RNA[, data := "Transcriptomic Only - ROC AUC: 0.74 ± 0.1"]
Prot[, data := "Proteomic Only - ROC AUC: 0.65 ± 0.19"]
WSI[, data := "Digital Pathology Only - ROC AUC: 0.57 ± 0.18"]
Int[, data := "Integrated - ROC AUC: 0.84 ± 0.1"]


All <- rbind(Clin,DNA,RNA,Prot,WSI,Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Clinic Only - ROC AUC: 0.61 ± 0.16","Genomic Only - ROC AUC: 0.79 ± 0.15",
                                    "Transcriptomic Only - ROC AUC: 0.74 ± 0.1","Proteomic Only - ROC AUC: 0.65 ± 0.19",
                                    "Digital Pathology Only - ROC AUC: 0.57 ± 0.18","Integrated - ROC AUC: 0.84 ± 0.1"))

ROC_TDM1 <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +
  scale_colour_manual(values = c("#434279", "#f2a104","#72a2c0","#9467bd","#00743f","#d62728")) +
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_TDM1

ggsave(ROC_TDM1, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/TDM1_roc_curves.pdf", width = 6.2, height = 6.2)
##########DHP##########
library(data.table);library(ggplot2);library(Polychrome)
library(data.table);library(ggplot2);library(Polychrome)
Int=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_0_Integrate/all_no_norm_TREAT_0_selfeats_10_run4_ROC_values.csv")
Clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_0_clinical/all_no_norm_TREAT_0_selfeats_6_clinical_ROC_values.csv")
DNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_0_DNA/all_no_norm_TREAT_0_selfeats_10_DNA_ROC_values.csv")
RNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_0_RNA/all_no_norm_TREAT_0_selfeats_10_RNA_ROC_values.csv")
Prot=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_0_Proteomics/all_no_norm_TREAT_0_selfeats_10_Proteomics_ROC_values.csv")
WSI=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_0_Image/all_no_norm_TREAT_0_selfeats_3_Image_ROC_values.csv")

rename <- c("fpr","tpr")

names(Int)[2:3] <- rename
names(Clin)[2:3] <- rename
names(DNA)[2:3] <- rename
names(RNA)[2:3] <- rename
names(Prot)[2:3] <- rename
names(WSI)[2:3] <- rename

Clin[, data := "Clinic Only - ROC AUC: 0.61 ± 0.15"]
DNA[, data := "Genomic Only - ROC AUC: 0.78 ± 0.14"]
RNA[, data := "Transcriptomic Only - ROC AUC: 0.84 ± 0.11"]
Prot[, data := "Proteomic Only - ROC AUC: 0.83 ± 0.11"]
WSI[, data := "Digital Pathology Only - ROC AUC: 0.65 ± 0.16"]
Int[, data := "Integrated - ROC AUC: 0.87 ± 0.1"]


All <- rbind(Clin,DNA,RNA,Prot,WSI,Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Clinic Only - ROC AUC: 0.61 ± 0.15","Genomic Only - ROC AUC: 0.78 ± 0.14",
                                    "Transcriptomic Only - ROC AUC: 0.84 ± 0.11","Proteomic Only - ROC AUC: 0.83 ± 0.11",
                                    "Digital Pathology Only - ROC AUC: 0.65 ± 0.16","Integrated - ROC AUC: 0.87 ± 0.1"))

ROC_DHP <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +
  scale_colour_manual(values = c("#434279", "#f2a104","#72a2c0","#9467bd","#00743f","#d62728"))+ 
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_DHP 
ggsave(ROC_DHP , file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/DHP_roc_curves.pdf", width = 6.2, height = 6.2)









