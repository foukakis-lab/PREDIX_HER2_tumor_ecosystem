############################
############################
############################
library(SeuratDisk);library(Seurat);library(tidyverse)
object=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/cell_state_integration/Seurat_cell_state_curated_localBPcell_Oct2025.rds")
table(object$cell_type)
library(MAST)

markers_wilcox <- FindAllMarkers(
  object,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.15,
  test.use        = "wilcox"
)
res=markers_wilcox%>%filter(avg_log2FC>1,pct.1>0.3,pct.2<0.2)
write.csv(res,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure7/cell_type_DEGs.csv",quote = F,row.names =F)


library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
df=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
df=df[df$sampleID=="1807_BL",]
df$group=df$cell_type_major_curated
df$group[df$cell_state=="Basal_SC"]="Basal_SC"
df$group[df$cell_state=="Her2E_SC"]="Her2E_SC"
df$group[df$cell_state=="LumA_SC"]="LumA_SC"
df$group[df$cell_state=="LumB_SC"]="LumB_SC"
df$cell_id=sub("^[^_]+_[^_]+_", "", row.names(df)) 
data_list <- split(df, df$sampleID)
# 循环导出每个子数据集
for (id in names(data_list)) {
  write.csv(data_list[[id]], file = paste0("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/",id, ".csv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
}


############### Export meta data for  Niche ################
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
data$cell_name=row.names(data)
threshold <- quantile(data$ADC_trafficking1[data$cell_state%in%c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC")], 0.75, na.rm = TRUE)
data$cell_state[data$cell_state%in%c("Her2E_SC")&data$ADC_trafficking1>threshold]="ADC trafficking Her2E_SC"  
table(data$cell_state)
write.csv(data,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.csv",quote = F,row.names =F)


# Xenium representative image, export cell label
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
df=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
nrow(df)
df$group=df$cell_state
df$cell_id=sub("^[^_]+_[^_]+_", "", row.names(df)) 
data_list <- split(df, df$sampleID)
# 循环导出每个子数据集
for (id in names(data_list)) {
  write.csv(data_list[[id]], file = paste0("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/sample_cellstate_label/",id, ".csv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
}

df$group=df$cell_type
data_list <- split(df, df$sampleID)
# 循环导出每个子数据集
for (id in names(data_list)) {
  write.csv(data_list[[id]], file = paste0("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/sample_majorcell_label/",id, ".csv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
}


df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baseline_niche_curated.csv")
df$merged_niche=paste0("Niche",df$merged_niche)
df$group=df$merged_niche
data_list <- split(df, df$sample)
# 循环导出每个子数据集
for (id in names(data_list)) {
  write.csv(data_list[[id]], file = paste0("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/sample_niche/",id, ".csv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
}


############################
######### Fig1b-c #########
############################
library(Seurat);library(scCustomize)
library(SeuratObject);library(BPCells)
library(dplyr);library(ggplot2);library(patchwork)
library(data.table)
options(future.globals.maxSize = 50 * 1024^3)  
gc()
object=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/cell_state_integration/Seurat_cell_state_curated_localBPcell_Oct2025.rds")
nrow(object@meta.data)
Idents(object) <- factor(Idents(object), 
                       levels =c("Epithelial", "Tcell", "Bcell",
                                     "Myeloid","Mesenchymal","Endothelial",
                                     "Adipocytes"))
col=c("#D1352B","#3C77AF","#FCED82",
      "#90C2E7","#7DBFA7","#9B5B33","grey")
p1=DimPlot(object,reduction = "umap.harmony3",cols=col,
           alph=0.7,label = T)+NoLegend()
p1
ggsave(p1, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Fig6b.pdf", width=5, height=5)

all.markers <- FindAllMarkers(object = object)
table(Idents(object))
res=all.markers%>%filter(avg_log2FC>1,pct.1>0.3,pct.2<0.2)
write.csv(res,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure7/cell_type_DEGs.csv",quote = F,row.names =F)

cell_type=data.frame(cell_name=colnames(object),
                     cell_type=object$cell_type,
                     cell_state=object$cell_state)
write.csv(cell_type,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure7/cell_type_state.csv",
          quote = F,row.names =F)

object@meta.data$group=NULL 
object$group=object$cell_type
# bubble plot for canonical markers
object$group <- factor(object$group, 
                         levels =rev(c("Epithelial", "Tcell", "Bcell",
                                       "Myeloid","Mesenchymal","Endothelial",
                                       "Adipocytes")))
Idents(object)=object$group
markers <- list(
  "Epithelial"=c("ERBB2","EPCAM","KRT19","ERBB3","ESR1"),
  "Tcell"=c("CD3E","CD8A","IL7R","CCL5","GZMK"),
  "Bcell"=c("TENT5C",'CD79A',"MZB1","FCRL5"),
  "Myeloid"=c("CD68","FCGR2A","MSR1","SLCO2B1"),
  "Mesenchymal"=c("ACTA2","COL11A1","PDGFRB","THBS2"),
  "Endothelial"=c("VWF","PLVAP","ENG"),
  "Adipocytes"=c("FABP4","PLIN1","PLIN4","ADIPOQ")
)

DotPlot(object,markers)+theme(axis.text.x = element_text(angle = 90))
# 10X4

### fig6d ###
data=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.csv")
data$patientID=as.character(data$patientID)
pid=fread("E:/Projects/PREDIX_HER2/Multimodal/Data_repository/PREDIX_HER2_multiomics_meta.txt")
pid$patientID=NA
pid$patientID=substr(pid$sampleID.wes,9,12)%>%as.character()
pid$patientID[is.na(pid$patientID)]=substr(pid$sampleID.rna[is.na(pid$patientID)],9,12)%>%as.character()
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
data=left_join(data,meta,by="patientID")%>%left_join(pid,by="patientID")
data$cell_type <- factor(data$cell_type, 
                         levels =c("Epithelial", "Tcell", "Bcell",
                                   "Myeloid","Mesenchymal","Endothelial",
                                   "Adipocytes"))
col=c("#D1352B","#3C77AF","#FCED82",
      "#90C2E7","#7DBFA7","#9B5B33","grey")
data=data[order(data$Arm),]
table(is.na(data$study_id))
table(data$study_id[data$Arm=="T-DM1"])
table(data$study_id[data$Arm=="DHP"])
data$study_id=factor(data$study_id,levels = c("P103","P115","P119","P133","P135",
                                              "P137","P150","P153","P186","P188","P189","P3","P33","P48","P5","P69","P75","P98",
                                              "P114","P132","P182","P184","P41","P43","P63","P84"))

mytable=data %>%
  group_by(study_id,cell_type) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
write.table(mytable,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Fig6d.txt",quote = F,row.names =F,sep="\t")

library(ggpubr)
ggbarplot(mytable, "study_id", "freq",
             fill = "cell_type", 
             color = "cell_type", 
             palette = col,
             legend = "right") +
  labs(fill = "cell type", color = "cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7))
#5.5X3


####################################
################ pCR ###############
####################################
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
table(data$cell_state)
tumor=data[data$cell_state%in%c("Her2E_SC","LumA_SC","LumB_SC","Basal_SC"),]
threshold <- quantile(tumor$ADC_trafficking1, 0.75, na.rm = TRUE)
tumor$cell_state[tumor$cell_state=="Her2E_SC"&tumor$ADC_trafficking1>threshold]="ADC trafficking Her2E_SC"  
table(tumor$cell_state)
tumor_cell_proportion <- tumor %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
library(tidyverse) 
tumor_cell_proportion_wide <- tumor_cell_proportion %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
tumor_cell_proportion_wide=left_join(tumor_cell_proportion_wide,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=tumor_cell_proportion_wide%>%select(c("Arm","Response","Her2E_SC","ADC trafficking Her2E_SC")) #"Basal_SC",,"LumA_SC","LumB_SC"
write.table(d,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Fig6e.txt",quote = F,row.names =F,sep="\t")
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Fig3c <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.99)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Fig3c

ggsave(Fig3c, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Fig6e.pdf", width=4, height=4)
######Immune and stromal####
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
immune=c("Bcell","Myeloid","Tcell")
Immune=data[data$cell_type%in%immune,]
table(Immune$cell_type,Immune$sampleID)
table(Immune$sampleID)
#Immune <- Immune[!Immune$sampleID %in% c("1228_BO"), ] #"1404_BL","1512_BL"
#Immune <- Immune[!Immune$sampleID %in% c("1404_BL"), ] # can include
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
d=Immune_wide%>%select(c("Arm","Response","APOC1_macrophages","CD8_Trm","CTLs","CXCR4_T","Cycling T","ECM_macrophages","LAM","Langerhans",          
    "M1_macrophages","M2_macrophages","MARCO_macrophages","Mast cells","Memory B","NK","NKT","Naive B","Tfh","Th17","Treg","cDC1","cDC2")) 
d=Immune_wide%>%select(c("Arm","Response","Mast cells","Treg")) 
# unique(data$cell_state[data$cell_type_major_curated%in%immune])
write.table(d,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Fig6f.txt",quote = F,row.names =F,sep="\t")
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Fig3c <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.99)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Fig3c
ggsave(Fig3c, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Fig6f.pdf", width=4, height=4)

###### stromal ######
immune=c("Endothelial","Mesenchymal")
Immune=data[data$cell_type%in%immune,]
table(Immune$cell_type,Immune$sampleID)
table(Immune$sampleID)
#Immune <- Immune[!Immune$sampleID %in% c("1228_BO"), ] #"1404_BL","1512_BL"
#Immune <- Immune[!Immune$sampleID %in% c("1404_BL"), ] # can include
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
d=Immune_wide%>%select(c("Arm","Response","Cycling","ECM-CAF","Inflammatory","Lymphatic",
                         "PVLs","Vascular","iCAF","myCAF")) 
# unique(data$cell_state[data$cell_type_major_curated%in%immune])
write.table(d,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Fig6f.txt",quote = F,row.names =F,sep="\t")
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Fig3c <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.99)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Fig3c

#############################################################
#######################Niche analysis########################
#############################################################
library(data.table);library(tidyverse);library(ggplot2)
options(future.globals.maxSize = 50 * 1024^3)  
gc()
data=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baseline_niche_curated.csv")  ## ? filtered cell
data$ID=paste0(data$sample,"_",data$cell_id)
df=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
df$ID=row.names(df)
df=left_join(df,data,by="ID")
table(is.na(df$cell_state))
table(df$merged_niche)
df$merged_niche=paste0("Niche",df$merged_niche)
table(df$cell_type_major_curated,df$merged_niche)
table(df$cell_state)
df$cell_type_major_curated[df$cell_state%in%c("LumA_SC","LumB_SC",
                                              "Her2E_SC","Basal_SC")]="Tumor"
table(df$merged_niche,df$sampleID)
table(df$sampleID)

df$cell_state=factor(df$cell_state,levels = c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC","Luminal progenitor","Myoepithelial",
                                              "Cycling T","Exhausted T","CTLs","CD8 Tem","Th1","Th2","Th17","Treg","NK",
                                              "Naive B","Memory B","Plasma","Cycling myeloid","Monocyte","cDC","pDC","LAM","M1 macrophages","M2 macrophages","MDSC","Neutrophils","Mast cells",
                                              "Cycling Endothelial","Lymphatic Endothelial","Vascular Endothelial","EMT-like CAF","iCAF","myCAF","PVLs","Adipocytes"))
table(is.na(df$cell_state))
df=df[!df$sampleID%in%c("1228_BO","1404_BL"),]
#####################################
#############recurrent niche#########
#####################################
niche_table <- table(df$merged_niche, df$sampleID)
niche_df <- as.data.frame.matrix(niche_table)
# 计算每个 niche 出现的样本数量
niche_df$Sample_Count <- rowSums(niche_df > 0)
# 计算总样本数量
total_samples <- ncol(niche_df) - 1  # 减去Sample_Count
# 计算出现频率
niche_df$Frequency <- niche_df$Sample_Count / total_samples
# 计算出现样本中的平均丰度
niche_df$Mean_Abundance <- apply(niche_df[, 1:total_samples], 1, function(x) mean(x[x > 0]))
# 定义 recurrent niche 规则
freq_threshold <- 0.8  # 30% 样本
abundance_threshold <- 1000  # 平均丰度大于100
niche_df$Recurrent <- (niche_df$Frequency >= freq_threshold) & (niche_df$Mean_Abundance >= abundance_threshold)
# 查看 recurrent niches
recurrent_niches <- niche_df[niche_df$Recurrent, ]
table(recurrent_niches$Recurrent)
row.names(recurrent_niches)
########define niche######
df=df[df$merged_niche%in%row.names(recurrent_niches),]
# 统计recurrent niche中，不同细胞类型的数量
niche_cell_counts <- df %>%
  group_by(merged_niche, cell_type_major_curated) %>%
  summarise(cell_count = n()) %>%
  tidyr::pivot_wider(names_from = cell_type_major_curated, values_from = cell_count, values_fill = 0)

niche_composition <- niche_cell_counts %>%
  mutate(Total = Adipocytes + `B cells` + `Endothelial cells` + `Epithelial cells` + `Mast cells` +
           `Mesenchymal cells` + `Myeloid cells` + Plasma + `T cells` + Tumor)

# 计算各类细胞占比
niche_composition <- niche_composition %>%
  mutate(
    Tumor_ratio = Tumor / Total,
    Immune_ratio = (`B cells` + `Myeloid cells` + Plasma + `T cells` + `Mast cells`) / Total,
    Stromal_ratio = (Adipocytes + `Endothelial cells` + `Mesenchymal cells`) / Total,
    Immune_Stromal_ratio = Immune_ratio + Stromal_ratio
  )


niche_composition <- niche_composition %>%
  mutate(Niche_Type = case_when(
    Tumor_ratio >= 0.3 & Immune_Stromal_ratio>0.2~ "Tumor-Immune/Stromal niche",   #& Immune_Stromal_ratio >= 0.3 
    Tumor_ratio < 0.3 ~ "Immune/Stromal niche",
    TRUE ~ "Other niche"
  ))
table(niche_composition$Niche_Type)
table(niche_composition$merged_niche,niche_composition$Niche_Type)
df$merged_niche=factor(df$merged_niche,levels = rev(c("Niche0","Niche1","Niche2","Niche3","Niche4",
                                                      "Niche10","Niche11","Niche14","Niche15",
                                                      "Niche5","Niche6","Niche7","Niche8","Niche12","Niche13","Niche16","Niche24","Niche41")) ) # Immune/Stromal niche
########cell proportion by Niche#######
library(viridis)
cell_proportion_by_Niche=df %>%
  group_by(cell_state, merged_niche) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cell_state) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
write.table(cell_proportion_by_Niche,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Fig6g.txt",quote = F,row.names =F,sep="\t")

cell_proportion_by_Niche$proportion[cell_proportion_by_Niche$proportion<0.05]=0
# 复制一列 proportion 作为绘图用
cell_proportion_by_Niche$plot_proportion <- cell_proportion_by_Niche$proportion
# 0的值替换为很小值，保证ggplot可以绘制
# 创建颜色分组（0和非0的颜色区分开）
cell_proportion_by_Niche$plot_group <- ifelse(cell_proportion_by_Niche$proportion == 0, "zero", "nonzero")

g1 = ggplot(cell_proportion_by_Niche, aes(x = cell_state, y = merged_niche, size = plot_proportion, color = plot_proportion)) +
  geom_point(data = subset(cell_proportion_by_Niche, plot_group == "zero"), color = "white", size = 1.5) +  # proportion = 0 的灰色点
  geom_point(data = subset(cell_proportion_by_Niche, plot_group == "nonzero"), aes(color = plot_proportion)) + # 非0用viridis
  scale_size(range = c(1, 10)) +
  scale_color_viridis(option = "D", direction = 1, name = "Proportion") + # viridis调色板
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Niche vs Cell Proportion", x = "Cell Type", y = "Niche", size = "Proportion")
g1

ggsave(g1, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Fig6g.pdf", width=11, height=6.5)


# cell state proportion across niche 
# 自动生成足够颜色
my_colors <- c("#1f78b4","#a6cee3","#fb9a99","#e31a1c", "#2CA02C", "#00FFC0", "#FF007C", "#1F77B4",
               "#E377C2", "#B800FF", "#8C564B", "#98DF8A", "#9467BD", "#D0FF00", "#C49C94", "#17BECF",
               "#C5B0D5", "#BCBD22", "#FF6400", "#6CFF00", "#9EDAE5", "#F7B6D2", "#C7C7C7", "#FF9896",
               "#00D8FF", "#000FFF", "#07FF00", "#FF00E0", "#5400FF", "#FFC800", "#DBDB8D", "#0074FF",
               "#FFBB78", "#D62728", "#7F7F7F", "#00FF5C")


df$merged_niche=factor(df$merged_niche,levels = c("Niche0","Niche1","Niche2","Niche3","Niche4",
                                                  "Niche10","Niche11","Niche14","Niche15",
                                                  "Niche5","Niche6","Niche7","Niche8","Niche12",
                                                  "Niche13","Niche16","Niche24","Niche41") ) # Immune/Stromal niche
mytable=df %>%
  group_by(merged_niche,cell_state) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
library(ggpubr)
p2=ggbarplot(mytable, "merged_niche", "freq",
             fill = "cell_state", 
             color = "cell_state", 
             palette =my_colors,
             legend = "right") +
  labs(fill = "cell_state", color = "cell_state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7))
p2
ggsave(p2, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/Figs10l.pdf", width=8, height=5)

####### Tumor-Immune/stromal Niche proportion ######### 
selected_niche=niche_composition$merged_niche[niche_composition$Niche_Type=="Tumor-Immune/Stromal niche"] 
niche_proportion <- df[df$merged_niche%in%selected_niche,]%>%
  group_by(sampleID, merged_niche) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sampleID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
niche_proportion_wide <- niche_proportion %>%
  select(sampleID, merged_niche, proportion) %>%
  pivot_wider(names_from = merged_niche, values_from = proportion, values_fill = 0)
niche_proportion_wide$patientID=substr(niche_proportion_wide$sampleID,1,4)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
niche_proportion_wide=left_join(niche_proportion_wide,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=niche_proportion_wide%>%select(c("Arm","Response"),matches("Niche"))
d1 <- reshape2::melt(d,id.vars=c("Arm","Response"))
####### Immune/stromal-stromal Niche proportion ######### 
selected_niche=niche_composition$merged_niche[niche_composition$Niche_Type=="Immune/Stromal niche"] 
niche_proportion <- df[df$merged_niche%in%selected_niche,]%>%
  group_by(sampleID, merged_niche) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sampleID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
niche_proportion_wide <- niche_proportion %>%
  select(sampleID, merged_niche, proportion) %>%
  pivot_wider(names_from = merged_niche, values_from = proportion, values_fill = 0)
niche_proportion_wide$patientID=substr(niche_proportion_wide$sampleID,1,4)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
niche_proportion_wide=left_join(niche_proportion_wide,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=niche_proportion_wide%>%select(c("Arm","Response"),matches("Niche"))
d2 <- reshape2::melt(d,id.vars=c("Arm","Response"))

d=rbind(d1,d2)
d0=d[!d$variable%in%c("Niche2","Niche6","Niche12"),]
d0$variable=factor(d0$variable,levels = c("Niche0","Niche1","Niche3","Niche4",
                                          "Niche10","Niche11","Niche14","Niche15",
                                          "Niche5","Niche7","Niche8",
                                          "Niche13","Niche16","Niche24","Niche41"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Fig3c <-
  ggplot(d0,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=2)+
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
Fig3c

ggsave(Fig3c, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/Figs10m.pdf", width=12, height=7)

# Niche12,2,6
d0=d[d$variable%in%c("Niche2","Niche6","Niche12"),]
write.table(d0,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Fig6i.txt",quote = F,row.names =F,sep="\t")
d0$variable=factor(d0$variable,levels = c("Niche2","Niche6","Niche12"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Fig3c <-
  ggplot(d0,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.95)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Fig3c

ggsave(Fig3c, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Fig6i.pdf", width=8, height=4)



Fig3c <-
  ggplot(d0,aes(x=Response,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.95)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Fig3c

ggsave(Fig3c, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Fig6i2.pdf", width=8, height=4)



