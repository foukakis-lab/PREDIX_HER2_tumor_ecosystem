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
table(is.na(data$cell_state))
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


df=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/Niche_meta.rds")
df$group=df$niche_clusters
table(is.na(df$cell_state))
data_list <- split(df, df$sampleID)
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

##########################################
######### Per cell distance #############
#########################################
library(Matrix)
library(ggplot2)
library(Seurat)
library(phenoptr)
library(dplyr)
library(purrr)
df=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/coordinates_curated.rds")
df$cell=df$cell_id
df=df[,c("cell","cell_state","cell_type")]
centroids=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/centroids_metadata.rds")
centroids
df=left_join(df,centroids)
df=df[!df$orig.ident%in%c("1310_BL","1807_BL","1204_BL"),] 

#df$cell_type[df$cell_state=="Mast cells"]="Mast_cell"
df$cell_type[df$cell_state%in%c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC")]="Tumor"
df=df[df$cell_type%in%c("Tumor","Bcell","Tcell","Myeloid"),]
row.names(df)=df$cell
df$sampleID=df$orig.ident
# For efficiently merging results

# 1. Get all unique sample IDs
samples <- unique(df$sampleID)

# 2. Create an empty list to temporarily store the results of each sample
dist_list <- list()

# 3. Loop through each sample to calculate distances
for (current_sample in samples) {
  
  message("Processing sample: ", current_sample)
  
  # Extract data for the current sample
  sample_df <- df %>% filter(sampleID == current_sample)
  
  # Check if the current sample is empty (in case it was completely filtered out earlier)
  if (nrow(sample_df) == 0) {
    message("  -> Current sample is empty, skipping.")
    next
  }
  
  # Convert to the specific format and column names strictly required by phenoptr
  cds <- sample_df %>%
    mutate(
      `Cell X Position` = x,
      `Cell Y Position` = y,
      Phenotype = as.character(cell_type),
      `Cell ID` = rownames(sample_df) 
    ) %>%
    select(`Cell X Position`, `Cell Y Position`, Phenotype, `Cell ID`) %>%
    as_tibble()
  
  # Calculate nearest distance directly with the prepared data (no need to filter again)
  dist_result <- tryCatch({
    find_nearest_distance(cds)
  }, error = function(e) {
    message("  -> Error calculating distance: ", e$message)
    return(NULL) # Return NULL if this sample causes an error
  })
  
  # If calculation is successful, record the sampleID and store it in the list
  if (!is.null(dist_result)) {
    dist_result$sampleID <- current_sample 
    dist_list[[current_sample]] <- dist_result
  }
}

# 4. Loop finished, combine the distance data from all samples into one data frame
final_dist_df <- bind_rows(dist_list)

saveRDS(final_dist_df,file='E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/ImmuneCell2tumor_distance.rds')


###########  K-distance plot #############
library(data.table);library(ggplot2)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
clin=clin[,c("patientID","Response","Arm")]
csd_with_distance=readRDS('E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/ImmuneCell2tumor_distance.rds')
df=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/coordinates_curated.rds")
df$patientID=NA
df$patientID=substr(df$cell_id,1,4)%>%as.integer()
df=left_join(df,clin,by="patientID")
#df$cell_type[df$cell_state=="Mast cells"]="Mast_cell"
df$cell_type[df$cell_state%in%c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC")]="Tumor"
df=df[df$cell_type%in%c("Tumor","Bcell","Tcell","Myeloid"),]
df=df[!df$patientID%in%c("1310","1807","1204"),] 
df=cbind(df,csd_with_distance)

d=df[df$cell_state%in%c("Treg"),]  # "Tcell","Myeloid"
ggplot(d[d$Arm=="T-DM1",], aes(`Distance to Tumor`, color=Response)) +
  geom_density(size=0.5)

