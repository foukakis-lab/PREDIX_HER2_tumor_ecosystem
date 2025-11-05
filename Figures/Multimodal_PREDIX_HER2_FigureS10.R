# FigureS10c
library(data.table);library(tidyverse);library(ggplot2)
options(future.globals.maxSize = 50 * 1024^3)  
gc()
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/Niche_meta.rds")
data=data[,c("cell.name","niche_clusters")]
df=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
df$cell_state[is.na(df$cell_state)]=df$cell_type[is.na(df$cell_state)]
table(df$cell_state)
df$ID=row.names(df)
df=left_join(df,data)
table(is.na(df$cell_state))
df$cell_state=factor(df$cell_state,levels = c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC","Normal epithelial",
                                              "Cycling T","CTLs","CD8_Trm","CXCR4_T","STAT1_T","Tfh","Th17","Treg","NK","NKT",
                                              "Memory B","CD24_B","Plasma","M1_macrophages","M2_macrophages","ECM_macrophages","APOC1_macrophages","MARCO_macrophages","LAM","Langerhans","cDC1","cDC2","pDC","Mast cells",
                                              "Cycling","Lymphatic","Vascular","Inflammatory","ECM-CAF","iCAF","myCAF","PVLs","Adipocytes"))
table(is.na(df$cell_state))
niche_table <- table(df$niche_clusters, df$sampleID)
########define niche######
niche_cell_counts <- df %>%
  group_by(niche_clusters, cell_type) %>%
  summarise(cell_count = n()) %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = cell_count, values_fill = 0)
########cell proportion by Niche#######
library(viridis)
cell_proportion_by_Niche=df %>%
  group_by(cell_state, niche_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cell_state) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

cell_proportion_by_patient=df %>%
  group_by(cell_state, patientID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cell_state) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

Niche_by_patient=df %>%
  group_by(niche_clusters, patientID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(niche_clusters) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()


write.table(Niche_by_patient,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/niche_clusters_by_patient.txt",quote = F,row.names =F,sep="\t")

cell_proportion_by_Niche$proportion[cell_proportion_by_Niche$proportion<0.05]=0
cell_proportion_by_Niche$plot_proportion <- cell_proportion_by_Niche$proportion
cell_proportion_by_Niche$plot_group <- ifelse(cell_proportion_by_Niche$proportion == 0, "zero", "nonzero")

cell_proportion_by_Niche$niche_clusters=factor(cell_proportion_by_Niche$niche_clusters,levels = rev(c("Niche1_Plasma","Niche2_Immunosuppressive",
      "Niche3_Stromal","Niche4_Immunoactive","Niche5_Luminal_tumor","Niche6_HER2E_tumor","Niche7_Basal_tumor","Niche8_Blood_vessel")))
g1 = ggplot(cell_proportion_by_Niche, aes(x = cell_state, y = niche_clusters, size = plot_proportion, color = plot_proportion)) +
  geom_point(data = subset(cell_proportion_by_Niche, plot_group == "zero"), color = "white", size = 1.5) +  # proportion = 0 的灰色点
  geom_point(data = subset(cell_proportion_by_Niche, plot_group == "nonzero"), aes(color = plot_proportion)) + # 非0用viridis
  scale_size(range = c(1, 10)) +
  scale_color_viridis(option = "D", direction = 1, name = "Proportion") + # viridis调色板
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Niche vs Cell Proportion", x = "Cell Type", y = "Niche", size = "Proportion")
g1

ggsave(g1, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS10/Figs10c.pdf", width=15, height=7.5)


# tumor Niche figS10d
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/Niche_meta.rds")
data=data[data$niche_clusters%in%c("Niche5_Luminal_tumor","Niche6_HER2E_tumor","Niche7_Basal_tumor"),]
data$patientID=substr(data$sampleID,1,4)
niche_proportion <- data%>%
  group_by(patientID, niche_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
niche_proportion_wide <- niche_proportion %>%
  select(patientID, niche_clusters, proportion) %>%
  pivot_wider(names_from = niche_clusters, values_from = proportion, values_fill = 0)
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
niche_proportion_wide=left_join(niche_proportion_wide,meta,by="patientID")
d=niche_proportion_wide%>%select(c("Arm","Response"),matches("Niche"))
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Fig <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.90)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Fig
ggsave(Fig, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS10/Figs10d.pdf", width=6, height=4)


# Immune/stromal figS10e
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/Niche_meta.rds")
data=data[data$niche_clusters%in%c("Niche1_Plasma","Niche2_iCAF","Niche3_Stromal","Niche4_Immunoactive","Niche8_Blood_vessel"),]
data$patientID=substr(data$sampleID,1,4)
niche_proportion <- data%>%
  group_by(patientID, niche_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
niche_proportion_wide <- niche_proportion %>%
  select(patientID, niche_clusters, proportion) %>%
  pivot_wider(names_from = niche_clusters, values_from = proportion, values_fill = 0)
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
niche_proportion_wide=left_join(niche_proportion_wide,meta,by="patientID")
d=niche_proportion_wide%>%select(c("Arm","Response"),matches("Niche"))
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Fig <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.90)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))

Fig

ggsave(Fig, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS10/Figs10e.pdf", width=9, height=4)










