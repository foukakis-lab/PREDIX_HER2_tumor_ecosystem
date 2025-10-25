# FigS10b-d 
##########Frequency of selected features##########
library(data.table);library(ggplot2);library(forcats);library(tidyverse)
All=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Both_Integrate/all_no_norm_selfeats_10_run4_Feature_weights.csv")
All=All[,-1];All$voting_percentage <- round(100*All$hybrid_ranks_ / sum(All$hybrid_ranks_),2)
All$group=sub("_.*", "",All$feat_name);All=All[1:20,]
DHP=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_0_Integrate/all_no_norm_TREAT_0_selfeats_10_run4_Feature_weights.csv")
DHP=DHP[,-1];DHP$voting_percentage <- round(100*DHP$hybrid_ranks_ / sum(DHP$hybrid_ranks_),2)
DHP$group=sub("_.*", "",DHP$feat_name);DHP=DHP[1:20,]
TDM1=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/Treat_1_Integrate/all_no_norm_TREAT_1_selfeats_10_run4_Feature_weights.csv")
TDM1=TDM1[,-1];TDM1$voting_percentage <- round(100*TDM1$hybrid_ranks_ / sum(TDM1$hybrid_ranks_),2)
TDM1$group=sub("_.*", "",TDM1$feat_name);TDM1=TDM1[1:20,]

group_colors <- c(
  Clin = "#434279", 
  DNA = "#f2a104", 
  RNA = "#72a2c0",
  Prot= "#9467bd",
  WSI = "#00743f"
)
#All
All$feat_name=factor(All$feat_name,levels = c("RNA_mRNA-PGR","WSI_Cell_Interaction","RNA_Fatty_acid_metabolism", 
                                              "RNA_mRNA-ESR1","Prot_CDK12","RNA_Glutathione_metabolism",
                                              "DNA_RAB11FIP1_CNA","RNA_sspbc.subtype_Her2","DNA_RPL19_CNA",             
                                              "Prot_SLC12A2","Prot_RPL19","Prot_ERBB2",                
                                              "RNA_Glycolysis","Prot_MIEN1","RNA_HER2DX_HER2_amplicon",  
                                              "Clin_ER","RNA_Treg","Prot_PPP1R1B",              
                                              "RNA_Exosome","Prot_EEA1")) 
All$feat_name <- fct_rev(All$feat_name)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
All_feature=ggplot(data=All,aes(feat_name,voting_percentage,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_manual(values = group_colors) +
  labs(x="",y="Voting percentage (%)")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"),
        legend.position = "top")+
  coord_flip()+
  scale_y_continuous(limits = c(0, 10)) # 设置X轴范围为0到100
All_feature
ggsave(All_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/All_feature_rank.pdf", width = 6, height = 5)

#DHP
DHP$feat_name = factor(DHP$feat_name, levels = c(
  "RNA_Th2 cells","Prot_RPL19","Prot_CDK12","Prot_MIEN1","Prot_ERBB2_PG",
  "RNA_Glutathione_metabolism","RNA_mRNA-PGR","Prot_ERBB2","RNA_mRNA-ERBB2","RNA_TAM_M2",
  "DNA_RPL19_CNA","Prot_EEA1","Prot_FLOT1","Prot_PPFIA1","RNA_Glycolysis","WSI_Cell_Interaction",
  "DNA_BRCA2_CNA","RNA_Taxane_response","RNA_FCGR3A","Prot_PPP1R1B"
))
DHP$feat_name <- fct_rev(DHP$feat_name)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
DHP_feature=ggplot(data=DHP,aes(feat_name,voting_percentage,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_manual(values = group_colors) +
  labs(x="",y="Voting percentage (%)")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"),
        legend.position = "top")+
  coord_flip()+
  scale_y_continuous(limits = c(0, 12)) 
DHP_feature
ggsave(DHP_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/DHP_feature_rank.pdf", width = 6, height = 5)

#T-DM1
TDM1$feat_name = factor(TDM1$feat_name, levels = c(
  "DNA_RAB11FIP1_CNA",
  "RNA_mRNA-PGR",
  "Prot_PPP1R1B",
  "Prot_SLC12A2",
  "DNA_BRCA2_CNA",
  "RNA_mRNA-ESR1",
  "DNA_CNV_burden",
  "RNA_Exosome",
  "Prot_EEA1",
  "RNA_Fatty_acid_metabolism",
  "DNA_coding_mutation_TP53_oncokb",
  "RNA_FCGR3B",
  "RNA_ABC_transporter",
  "Prot_ARL1",
  "RNA_MHC.I_19272155",
  "Prot_VAMP3",
  "RNA_Mast-cells",
  "WSI_Cell_Interaction",
  "RNA_HER2DX_luminal",
  "DNA_PIK3CA_CNA"
))
TDM1$feat_name <- fct_rev(TDM1$feat_name)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
TDM1_feature=ggplot(data=TDM1,aes(feat_name,voting_percentage,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_manual(values = group_colors) +
  labs(x="",y="Voting percentage (%)")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"),
        legend.position = "top")+
  coord_flip()+
  scale_y_continuous(limits = c(0, 8)) 
TDM1_feature
ggsave(TDM1_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/TDM1_feature_rank.pdf", width = 7, height = 5)

# FigS10e-f

library(data.table);library(ggplot2);library(Polychrome)
RNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/Best_ML_Results_12_10_25/BestResults/RNA_EXTERNAL_TREAT_0.csv")

rename <- c("fpr","tpr")
names(RNA)[2:3] <- rename
RNA[, data := "Transcriptomic Only - ROC AUC: 0.76"]

ROC_TDM1 <- ggplot(RNA, aes(y = tpr, x = fpr, colour = data, group = data)) + 
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




