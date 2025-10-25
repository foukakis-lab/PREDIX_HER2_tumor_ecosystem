# FigureS9a
library(Seurat);library(scCustomize)
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
data=data[data$cell_state%in%c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"),]
data <- data %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
data <- data %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
pid=fread("E:/Projects/PREDIX_HER2/Multimodal/Data_repository/PREDIX_HER2_multiomics_meta.txt")
pid$patientID=NA
pid$patientID=substr(pid$sampleID.wes,9,12)%>%as.character()
pid$patientID[is.na(pid$patientID)]=substr(pid$sampleID.rna[is.na(pid$patientID)],9,12)%>%as.character()
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
data=left_join(data,meta,by="patientID")%>%left_join(pid,by="patientID")%>%left_join(rna,by="patientID")
table(data$sspbc.subtype)
data$sspbc=NA
#data$sspbc[data$sspbc.subtype%in%c("LumA","LumB")]="LumAorB"
#data$sspbc[data$sspbc.subtype%in%c("Her2")]="Her2"

d=data%>%select(c("sspbc.subtype","LumA_SC","LumB_SC","Her2E_SC","Basal_SC")) #"Basal_SC",,"LumA_SC","LumB_SC"
d$sspbc.subtype=as.character(d$sspbc.subtype)
d$sspbc.subtype[d$sspbc.subtype%in%c("LumA","LumB")]="LumAorB"
d <- reshape2::melt(d,id.vars=c("sspbc.subtype"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
FigureS9a <-
  ggplot(d,aes(x=variable,y=value,fill=variable))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(sspbc.subtype~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_manual(values = c(
    "#1f78b4","#a6cee3","#fb9a99","#e31a1c"
  ))+
  stat_compare_means(aes(group=variable),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.95)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
FigureS9a

# FigureS9b Tumor cell
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
tumor=data[data$cell_state%in%c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC"),]
data$ADC_trafficking1
table(tumor$cell_state)
table(tumor$cell_state,tumor$sampleID)
tumor_cell_proportion <- tumor %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
tumor_cell_proportion_wide <- tumor_cell_proportion %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
tumor_cell_proportion_wide=left_join(tumor_cell_proportion_wide,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=tumor_cell_proportion_wide%>%select(c("Arm","Response","LumA_SC","LumB_SC","Her2E_SC","Basal_SC")) #
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
FigureS9b <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.98)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
FigureS9b
# 8X4

# FigureS9c Immune cell
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
data$cell_state=factor(data$cell_state,levels =c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC","Normal epithelial",
            "Cycling T","CTLs","CD8_Trm","CXCR4_T","STAT1_T","Tfh","Th17","Treg","NK","NKT",
              "Memory B","CD24_B","Plasma","M1_macrophages","M2_macrophages","ECM_macrophages","APOC1_macrophages","MARCO_macrophages","LAM","Langerhans","cDC1","cDC2","pDC","Mast cells",
              "Cycling","Lymphatic","Vascular","Inflammatory","ECM-CAF","iCAF","myCAF","PVLs") )

immune=c("Tcell","Myeloid","Bcell")
Immune=data[data$cell_type%in%immune,]
Immune_proportion <- Immune %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
Immune_wide <- Immune_proportion %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
Immune_wide =left_join(Immune_wide ,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=Immune_wide%>%select(c("Arm","Response","Cycling T","CTLs","CD8_Trm","CXCR4_T","STAT1_T","Tfh","Th17","Treg","NK","NKT",
                         "Memory B","CD24_B","Plasma","M1_macrophages","M2_macrophages","ECM_macrophages","APOC1_macrophages","MARCO_macrophages","LAM","Langerhans","cDC1","cDC2","pDC")) # 
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
str(d)
FigureS9c <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=3)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.98)+
  labs(y="Immune proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
FigureS9c
ggsave(FigureS9c, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS9/FigureS9c.pdf", width=20, height=8)

#FigureS9d
Endo=c("Endothelial")
Endo=data[data$cell_type%in%Endo,]
Endo_proportion <- Endo %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
Endo_wide <- Endo_proportion %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
Endo_wide =left_join(Endo_wide ,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=Endo_wide%>%select(c("Arm","Response","Cycling","Lymphatic","Vascular","Inflammatory")) # 
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
str(d)
FigureS9d <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.95)+
  labs(y="Endothelial cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
FigureS9d
ggsave(FigureS9d, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS9/FigureS9d.pdf", width=8, height=3)


#FigureS9d
stromal=c("Mesenchymal")
stromal=data[data$cell_type%in%stromal,]
stromal_proportion <- stromal %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
stromal_wide <- stromal_proportion %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
stromal_wide =left_join(stromal_wide ,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=stromal_wide%>%select(c("Arm","Response","ECM-CAF","iCAF","myCAF","PVLs")) # 
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
str(d)
FigureS9e <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.95)+
  labs(y="Mesenchymal stromal cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
FigureS9e
ggsave(FigureS9e, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS9/FigureS9e.pdf", width=8, height=3)
