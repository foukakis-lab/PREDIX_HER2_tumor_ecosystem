# fig7a
#########feature correlation##########
library(data.table);library(ggplot2);library(ggpubr);library(tidyverse);library(ComplexHeatmap)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
#data=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
data=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2_withNA.txt")%>%as.data.frame()
data$WSI_Immune_Cell_prop=as.numeric(data$WSI_Immune_Cell_prop)
data$WSI_Distance_tumor_immune=as.numeric(data$WSI_Distance_tumor_immune)
data$WSI_Cell_Interaction=as.numeric(data$WSI_Cell_Interaction)
library(fastDummies)
data <- dummy_cols(data, select_columns = "RNA_sspbc.subtype")
data$RNA_sspbc.subtype=NULL
data$RNA_sspbc.subtype_NA=NULL
data$Clin_prolifvalu=as.numeric(data$Clin_prolifvalu)
data <- data %>%
  mutate(Clin_TUMSIZE = case_when(
    Clin_TUMSIZE == "<=20"    ~ 1,
    Clin_TUMSIZE == "21-50"   ~ 2,
    Clin_TUMSIZE == ">50"     ~ 3,
    Clin_TUMSIZE == "Unknown" ~ NA
  ))
data <- data %>%
  mutate(Clin_ANYNODES = case_when(
    Clin_ANYNODES == "N0"    ~ 1,
    Clin_ANYNODES == "N+"   ~ 2
  ))

clin_cols <- c("Clin_prolifvalu","Clin_TUMSIZE","Clin_ANYNODES")
rna_cols <- grep("^RNA_", colnames(data), value = TRUE)
dna_cols <- grep("^DNA_", colnames(data), value = TRUE)
prot_cols <- grep("^Prot_", colnames(data), value = TRUE)
wsi_cols <- grep("^WSI_", colnames(data), value = TRUE)
data$Arm=data$Clin_Arm


# univariate association---Clin
df=data[,c(clin_cols,"pCR","Arm")]
df <- df[complete.cases(df), ]
results_clin=Logistic_batch_uni(df,"pCR","Arm",clin_cols)%>%as.data.frame()
results=results_clin
sig_vars <- data.frame(
  biomarker = results$biomarker,
  whole_sig = results$whole_lr_p < 0.05,
  DHP_sig   = results$DHP_lr_p < 0.05,
  TDM1_sig  = results$TDM1_lr_p < 0.05
) %>%
  filter(whole_sig | DHP_sig | TDM1_sig) %>%
  pull(biomarker)
results_clin=results[results$biomarker%in%sig_vars,]

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
fig7a=results
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
pdf("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Fig7a_new.pdf", width = 9, height = 8)
draw(p)  
dev.off()

data=as.data.frame(data)
tab_csv=data[,c("Clin_Arm","Clin_TUMSIZE","Clin_ANYNODES","Clin_ER","Clin_prolifvalu",colnames(mat))]
tab_csv=na.omit(tab_csv)
write.csv(tab_csv, file = "E:/Projects/PREDIX_HER2/Multimodal/Table/supplemental_tab8.csv")
# fig7b-c

library(data.table);library(ggplotify);library(ggpubr)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/consensus_signatures.csv") 
df=df[df$Scenario=="DHP",]
df$Modality=factor(df$Modality,levels = c("Clin","DNA","RNA","Prot","WSI"))
group_colors <- c(
  Clin = "#434279", 
  DNA = "#f2a104", 
  RNA = "#72a2c0",
  Prot= "#9467bd",
  WSI = "#00743f"
)
df=df[order(df$Modality,-df$mean_imp),]
df$Feature <- factor(df$Feature, levels = rev(unique(df$Feature)))

p <- ggplot(df, aes(x = mean_imp, y = Feature, fill = Modality)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, alpha = 0.9) +
  scale_fill_manual(values = group_colors) + 
  labs(x = "Mean Importance", y = "") +     
  theme_classic() + 
  theme(
    axis.text.x = element_text(color = "black"), 
    axis.text.y = element_text(color = "black"), 
    legend.position = "top", 
    legend.title = element_text(face = "bold")
  )
p
Feature=df$Feature

library(IOBR);library(data.table);library(tableone)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
df$Clin_prolifvalu[df$Clin_prolifvalu=="Unknown"]=NA
df$Clin_prolifvalu=as.numeric(df$Clin_prolifvalu)
df$DNA_HRD=as.numeric(df$DNA_HRD)
df=df[df$Clin_Arm=="T-DM1",]
Feature=setdiff(Feature,c("Clin_ER","Clin_ANYNODES","Clin_TUMSIZE","Prot_ERBB2_PG","RNA_sspbc_LumA"))
str(df[,Feature]) 
res=batch_wilcoxon(
  df,
  target = "pCR",
  feature =c('DNA_ERBB2_CNA',"DNA_BRCA2_CNA")
)
res=res[res$p.value<0.05,]

df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
df=df[df$Clin_Arm=="DHP",]
cat=c("Clin_ER","Clin_ANYNODES","Clin_TUMSIZE","Prot_ERBB2_PG","RNA_sspbc_LumA")
tableOne <- CreateTableOne(vars = cat, strata = c("pCR"), data = df,
                           factorVars = cat)
print(tableOne)


df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/consensus_signatures.csv") 
df=df[df$Scenario=="T-DM1",]
df$Modality=factor(df$Modality,levels = c("Clin","DNA","RNA","Prot","WSI"))
group_colors <- c(
  Clin = "#434279", 
  DNA = "#f2a104", 
  RNA = "#72a2c0",
  Prot= "#9467bd",
  WSI = "#00743f"
)
df=df[order(df$Modality,-df$mean_imp),]
df$Feature <- factor(df$Feature, levels = rev(unique(df$Feature)))

p <- ggplot(df, aes(x = mean_imp, y = Feature, fill = Modality)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, alpha = 0.9) +
  scale_fill_manual(values = group_colors) + 
  labs(x = "Mean Importance", y = "") +     
  theme_classic() + 
  theme(
    axis.text.x = element_text(color = "black"), 
    axis.text.y = element_text(color = "black"), 
    legend.position = "top", 
    legend.title = element_text(face = "bold")
  )

Feature=df$Feature
print(p)
ggsave("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Fig7c.pdf", plot = p, width = 6, height = 6, dpi = 300)

library(IOBR);library(data.table);library(tableone)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
df$Clin_prolifvalu[df$Clin_prolifvalu=="Unknown"]=NA
df$Clin_prolifvalu=as.numeric(df$Clin_prolifvalu)
df$DNA_HRD=as.numeric(df$DNA_HRD)
df=df[df$Clin_Arm=="T-DM1",]
Feature=setdiff(Feature,c("Clin_ER","Clin_ANYNODES","Clin_TUMSIZE","Prot_ERBB2_PG","RNA_sspbc_LumB"))
str(df[,Feature]) 
res=batch_wilcoxon(
  df,
  target = "pCR",
  feature = Feature
)
res=res[res$p.value<0.05,]

df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
df=df[df$Clin_Arm=="T-DM1",]
cat=c("Clin_ER","Clin_ANYNODES","Clin_TUMSIZE","Prot_ERBB2_PG","RNA_sspbc_LumA")
tableOne <- CreateTableOne(vars = cat, strata = c("pCR"), data = df,
                           factorVars = cat)
print(tableOne)


# fig7d-e
library(data.table);library(ggplot2);library(ggpubr);library(tidyverse);library(ComplexHeatmap)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
#data=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
data=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2_withNA.txt")%>%as.data.frame()
data=data[!is.na(data$`RNA_mRNA-ERBB2`),]
data=data[!is.na(data$DNA_COSMIC.Signature.6),]

library(dplyr)

# Step-down stratification for S1, S2, and S3 in a single step
data <- data %>%
  mutate(
    Group = case_when(
      # Step 1: Identify Group S1 (Benefit from traditional dual blockade + chemotherapy)
      # Check if variables are >= their respective medians
      (RNA_Taxane_response >= median(RNA_Taxane_response, na.rm = TRUE) | 
        DNA_LOH_Del_burden >= median(DNA_LOH_Del_burden, na.rm = TRUE)) & 
        DNA_ERBB2_CNA >= quantile(DNA_ERBB2_CNA, probs = 0.75, na.rm = TRUE) #median(DNA_ERBB2_CNA, na.rm = TRUE)
      ~ "S1",
      
      # Step 2: From the remaining patients, identify Group S2 (Benefit from ADC)
      # Note: case_when evaluates sequentially. We don't need 'Group == "Other"' 
      # because it only evaluates patients not already assigned to S1.
      (RNA_ADC_trafficking >= median(RNA_ADC_trafficking, na.rm = TRUE)) & 
        (`RNA_Mast-cells` <= median(`RNA_Mast-cells`, na.rm = TRUE)| RNA_Neutrophils <= median(RNA_Neutrophils, na.rm = TRUE))
      ~ "S3",
      
      # Step 3: All remaining patients who do not meet the above criteria are assigned to Group S3
      TRUE ~ "S2"
    )
  )

# Check the final grouping distribution
print("Total cohort grouping statistics:")
table(data$Group)

# Cross-tabulation check with clinical outcome (pCR) and treatment arm
print("Cross-tabulation of Group vs pCR vs Treatment Arm:")
table(data$Group, data$pCR, data$Clin_Arm)


whole<- glm(as.numeric(pCR) ~ Clin_Arm+Clin_ER, family = "binomial", data = data[data$Group=="S1",])
ShowRegTable(whole)
whole<- glm(as.numeric(pCR) ~ Clin_Arm+Clin_ER, family = "binomial", data = data[data$Group=="S2",])
ShowRegTable(whole)
whole<- glm(as.numeric(pCR) ~ Clin_Arm+Clin_ER, family = "binomial", data = data[data$Group=="S3",])
ShowRegTable(whole)

interaction_1<- glm(as.numeric(pCR) ~ Clin_Arm+Clin_ER+Group, family = "binomial", data = data)
interaction_2<- glm(as.numeric(pCR) ~ Clin_Arm+Clin_ER+Group+Group*Clin_Arm, family = "binomial", data = data)
lrtest(interaction_1,interaction_2)


library(tidyverse);library(tableone);library(data.table);library(lmtest);library(forestploter)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Subgroup_Forestplot.csv")
fig7f=df
df$DHP <- ifelse(is.na(df$DHP), "", df$DHP)
df$`T-DM1` <- ifelse(is.na(df$`T-DM1`), "", df$`T-DM1`)
df$`P for interaction` <- ifelse(is.na(df$`P for interaction`), "", df$`P for interaction`)
df$` ` <- paste(rep(" ", 10), collapse = " ")
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
            ticks_at = c(0,0.5, 1, 2, 3),
            theme = tm)

# Print plot
plot(p)



library(openxlsx)
data_list=list("fig7a"=fig7a%>%as.data.frame(),
               "fig7b"=fig7b%>%as.data.frame(),
               "fig7c"=fig7c%>%as.data.frame(),
               "fig7d"=fig7d%>%as.data.frame())
openxlsx::write.xlsx(data_list,file='E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Figure7.xlsx')




