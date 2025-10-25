##############
# FigureS8a,b #
##############
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
pid=fread("E:/Projects/PREDIX_HER2/Multimodal/Data_repository/PREDIX_HER2_multiomics_meta.txt")
pid$patientID=NA
pid$patientID=substr(pid$sampleID.wes,9,12)%>%as.character()
pid$patientID[is.na(pid$patientID)]=substr(pid$sampleID.rna[is.na(pid$patientID)],9,12)%>%as.character()
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
data=left_join(data,meta,by="patientID")%>%left_join(pid,by="patientID")
data=data[order(data$Arm),]
table(is.na(data$study_id))
table(data$study_id[data$Arm=="T-DM1"])
table(data$study_id[data$Arm=="DHP"])

meta_summary <- data %>%
  group_by(study_id,patientID) %>%
  summarise(num_cells_detected = n(), .groups = "drop")%>%left_join(meta,by="patientID")
meta_summary$study_id=factor(meta_summary$study_id,levels = c("P103","P115","P119","P133","P135",
                                                              "P137","P150","P153","P186","P188","P189","P3","P33","P48","P5","P69","P75","P98",
                                                              "P114","P132","P182","P184","P41","P43","P63","P84"))

p1 <- ggplot(meta_summary, aes(x = study_id, y = num_cells_detected)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#3C77AF") +
  labs(title = "Number of cells detected", x = "sample ID", y = "cells detected") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  annotate("text",x="P133", y = 120000, 
           label = paste("Median=",median(meta_summary$num_cells_detected)), 
           size = 4, color = "black")
p1

data=fread("E:/Projects/PREDIX_HER2_longitudinal/data/Xenium/merged_metrics_summary.csv")
data$patientID=substr(data$region_name,1,4)
data$tpt=substr(data$region_name,6,7)
data$tpt[data$tpt%in%c("BL","BO")]="pre"
data=data[data$tpt=="pre",]
pid=fread("E:/Projects/PREDIX_HER2/Multimodal/Data_repository/PREDIX_HER2_multiomics_meta.txt")
pid$patientID=NA
pid$patientID=substr(pid$sampleID.wes,9,12)%>%as.character()
pid$patientID[is.na(pid$patientID)]=substr(pid$sampleID.rna[is.na(pid$patientID)],9,12)%>%as.character()
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
data=left_join(data,meta,by="patientID")%>%left_join(pid,by="patientID")
data=data[data$study_id!="P108",]
data$study_id=factor(data$study_id,levels = c("P103","P115","P119","P133","P135",
                                              "P137","P150","P153","P186","P188","P189","P3","P33","P48","P5","P69","P75","P98",
                                              "P114","P132","P182","P184","P41","P43","P63","P84"))
p2 <- ggplot(data, aes(x = study_id, y = median_transcripts_per_cell)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#3C77AF") +
  labs(title = "Median number of transcripts detected per cell", x = "sample ID", y = "Median number of transcripts per Cell") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  annotate("text",x="P133", y = 500, 
           label = paste("Median (total)=",median(data$median_transcripts_per_cell)), 
           size = 4, color = "black")
p2

p=p1+p2
p
ggsave(p, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS8/FigS8a_b.pdf", width=15, height=5)
# 12X5

################
# FigureS8c-h #
################
library(Seurat);library(scCustomize)
options(future.globals.maxSize = 50 * 1024^3)  
object=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/cell_state_integration/Seurat_cell_state_curated_localBPcell_Oct2025.rds")
object$cell_state[object$cell_state%in%c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC")]="Tumor epithelial"
object=subset(object,idents = c("Adipocytes"), invert = TRUE)
Idents(object)=object$cell_state
table(Idents(object))
Idents(object)=factor(Idents(object),levels =c("Tumor epithelial","Normal epithelial",
            "Cycling T","CTLs","CD8_Trm","CXCR4_T","STAT1_T","Tfh","Th17","Treg","NK","NKT",
            "Memory B","CD24_B","Plasma",
            "M1_macrophages","M2_macrophages","ECM_macrophages","APOC1_macrophages","MARCO_macrophages","LAM","Langerhans","cDC1","cDC2","pDC","Mast cells",
            "Cycling","Lymphatic","Vascular","Inflammatory",
            "ECM-CAF","iCAF","myCAF","PVLs") )
Epi=subset(object,idents = c("Tumor epithelial","Normal epithelial"))
Tcell=subset(object,idents = c("Cycling T","CTLs","CD8_Trm","CXCR4_T","STAT1_T","Tfh","Th17","Treg","NK","NKT"))
Bcell=subset(object,idents = c("Memory B","CD24_B","Plasma"))
Myeloid=subset(object,idents = c("M1_macrophages","M2_macrophages","ECM_macrophages","APOC1_macrophages","MARCO_macrophages","LAM","Langerhans","cDC1","cDC2","pDC","Mast cells"))
Endo=subset(object,idents = c("Cycling","Lymphatic","Vascular","Inflammatory"))
Mesenchymal=subset(object,idents = c("ECM-CAF","iCAF","myCAF","PVLs"))
table(object$cell_state[object$cell_type=="Tcell"])
library(patchwork)
# Epithelial
features = c("ACTA2","MKI67","EPCAM","TP63","KIT","KRT19","ERBB2","ESR1","PGR","FOXA1")
p=Clustered_DotPlot(seurat_object =Epi, features = features,
                     show_ident_legend = FALSE)
# 5X5
# Tcell
features = c(
  "TRAC","CD3E","CD8A","CD4","MKI67","TK1","TUBB","FCGR3A","KLRK1","NCAM1",
  "RORC","CD40LG","CCR6","KLRB1","DPP4","GZMA","GZMK","KLRK1","CCL4",
  "MIAT","ITGA1","ITGAE","FOXP3","SELL","CCR7","IL7R","CXCR4","PDCD1",
  "CXCL13","CXCR5","BCL6","ICOS","STAT1"
)
Clustered_DotPlot(seurat_object =Tcell, features = features,
                  show_ident_legend = FALSE)


# Bcell
features = c(
  "MS4A1","CXCR4","CD79A","XBP1","MZB1","TENT5C","TNFRSF17","CD24"
)
Clustered_DotPlot(seurat_object =Bcell, features = features,
                  show_ident_legend = FALSE)

# Myeloid cell
features = c(
  "CD163","IL1B","TNF","CCL4","CCL8","CXCL2","OSM","APOC1","MARCO",
  "CLEC9A","CD1C","CD1D","CLEC4C","IL3RA","CD207","CD93",
  "SPP1","CD36","LPL","PLIN2","CTSK","MMP9","ACP5","ITGB3","NFATC1",
  "CPA3","MS4A2","FCGR3A","CD14","CD68","CD86","CD163","TREM2","FABP4",
  "ITGAX","CLEC4C","CD33","FCGR3B","FPR2","CSF3R"
)
Clustered_DotPlot(seurat_object =Myeloid, features = features,
                  show_ident_legend = FALSE)

# Endothelial cell
features = c(
  "PECAM1", "VWF", "PROX1","MKI67","TK1","CXCL9","CXCL10","CXCL11"
)
Clustered_DotPlot(seurat_object =Endo, features = features,
                  show_ident_legend = FALSE)


# Mesenchymal cell
features=c("MCAM","PDGFRB","CXCL12","CXCL16","FN1","FAP","ACTA2","THBS1","CCN1","CCN2","CCN3")
Clustered_DotPlot(seurat_object = Mesenchymal, features = features,
                  show_ident_legend = FALSE)


