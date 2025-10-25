#FigureS7a
library(tidyverse)
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
data$manuTILS=as.numeric(data$manuTILS)
row.names(data)=data$patientID
data[,variable]=scale(data[,variable])
data=data[,c("manuTILS",variable)]
data=na.omit(data)
target_var <- data[, 1]  # 选择第一个变量作为目标变量
# 使用 sapply 函数计算 Spearman 相关系数和 p 值
cor_RNA <- sapply(data[, -1], function(x) {
  test <- cor.test(target_var, x, method = "spearman")
  c(correlation = test$estimate, p_value = test$p.value)
})%>%t()%>%as.data.frame()

#DNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
data$manuTILS=as.numeric(data$manuTILS)
variable=c("meanHED","TCRA.tcell.fraction.adj")
data[,variable]=scale(data[,variable])
data=data[,c("manuTILS",variable)]%>%as.data.frame()
data=na.omit(data)
target_var <- data[, 1]  # 选择第一个变量作为目标变量
# 使用 sapply 函数计算 Spearman 相关系数和 p 值
cor_DNA <- sapply(data[, -1], function(x) {
  test <- cor.test(target_var, x, method = "spearman")
  c(correlation = test$estimate, p_value = test$p.value)
})%>%t()%>%as.data.frame()

#integrative
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
data$manuTILS=as.numeric(data$manuTILS)
variable=c("Neoantigen")
data[,c("Neoantigen")]=scale(data[,c("Neoantigen")])
data=data[,c("manuTILS",variable)]%>%as.data.frame()
data=na.omit(data)
target_var <- data[, 1]  # 选择第一个变量作为目标变量
# 使用 sapply 函数计算 Spearman 相关系数和 p 值
test <- cor.test(target_var, data$Neoantigen, method = "spearman")
cor_integrated=data.frame(correlation.rho = test$estimate, p_value = test$p.value)
row.names(cor_integrated)="Neoantigen"
#image
library(data.table)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
image=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Digital_image/handcraft_features_FULL.csv")%>%as.data.frame()
image$V1=NULL
image$TREAT=NULL
image$ERPRdic=NULL
image$pCR=NULL
image=image[,c("patientID","CellProps__ALL__ImmuneCells","MinDist__ALL__Tumor__ImmuneCells__Mean","median_FeD_Cent_mst__edgelength_min_div_max__ALL")]
data=left_join(image,clin,by="patientID")%>%as.data.frame()
varibale=c("CellProps__ALL__ImmuneCells","MinDist__ALL__Tumor__ImmuneCells__Mean","median_FeD_Cent_mst__edgelength_min_div_max__ALL")
data=data[,c("manuTILS",varibale)]
data$manuTILS=as.numeric(data$manuTILS)
colnames(data)=c("manuTILS","Immune_Cell_prop","Distance_tumor_immune","Cell_Interaction")
variable=c("Immune_Cell_prop","Distance_tumor_immune","Cell_Interaction")
data[,variable]=scale(data[,variable])
data=data[,c("manuTILS",variable)]%>%as.data.frame()
data=na.omit(data)
target_var <- data[, 1]  # 选择第一个变量作为目标变量
# 使用 sapply 函数计算 Spearman 相关系数和 p 值
cor_image <- sapply(data[, -1], function(x) {
  test <- cor.test(target_var, x, method = "spearman")
  c(correlation = test$estimate, p_value = test$p.value)
})%>%t()%>%as.data.frame()
######merge#######
df=rbind(cor_RNA,cor_DNA)%>%rbind(cor_integrated)%>%rbind(cor_image)
df$correlation.rho=round(df$correlation.rho,digits = 2)
df=df[order(df$cor,decreasing = T),]
df$name=row.names(df)
df$p_value=format(df$p_value, digits = 3, scientific = TRUE)
df$p=paste0("P=",df$p_value)
colnames(df)=c("rho","p","name","p_value")
Immunogenicity=df[c("Neoantigen","meanHED","MHC.I_19272155"),]
Immune_Activation=df[c("B-cells","DC","Macrophages","T-cells","CD8-T-cells","Cytotoxic-cells","NK-cells","CD45","TILs","Immune_Cell_prop",
                       "TCR_clonality","BCR_clonality","TCRA.tcell.fraction.adj","FCGR3A","FCGR3B"),]
Immunosuppression=df[c("Dysfunction","Exclusion","CAF","Neutrophils","Mast-cells","Treg","TAM_M2","Th2 cells"),]
Spatial_TME=df[c("Distance_tumor_immune","Cell_Interaction"),]
library(ggplot2)
library(dplyr)
# 生成热图
p1=ggplot(Immunogenicity, aes(x = "", y = reorder(name, rho), fill = rho)) +
  geom_tile(color = "white", width = 1,height = 1) +  # 绘制热图
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),       # 指定颜色
    values = scales::rescale(c(-0.5, 0, 0.7)),    # -1映射为蓝色，0为白色，1为红色
    limits = c(-0.5, 0.7)                         # 设置颜色范围
  ) +
  labs(x = NULL, y = NULL, fill = "Correlation") +
  geom_text(aes(label = rho), color = "black",size=2.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  expand_limits(x = c(0,10)) +  # 扩展 x 轴的范围以显示 p_value
  geom_text(aes(x = 1.5, label = p_value), hjust = 0, nudge_x = 0.3, color = "black")+coord_fixed() # 热图右侧显示 p_value
p2=ggplot(Immune_Activation, aes(x = "", y = reorder(name, rho), fill = rho)) +
  geom_tile(color = "white", width = 1,height = 1) +  # 绘制热图
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),       # 指定颜色
    values = scales::rescale(c(-0.5, 0, 0.7)),    # -1映射为蓝色，0为白色，1为红色
    limits = c(-0.5, 0.7)                         # 设置颜色范围
  ) +
  labs(x = NULL, y = NULL, fill = "Correlation") +
  geom_text(aes(label = rho), color = "black",size=2.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  expand_limits(x = c(0,10)) +  # 扩展 x 轴的范围以显示 p_value
  geom_text(aes(x = 1.5, label = p_value), hjust = 0, nudge_x = 0.3, color = "black")+coord_fixed() # 热图右侧显示 p_value
p3=ggplot(Immunosuppression, aes(x = "", y = reorder(name, rho), fill = rho)) +
  geom_tile(color = "white", width = 1,height = 1) +  # 绘制热图
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),       # 指定颜色
    values = scales::rescale(c(-0.5, 0, 0.7)),    # -1映射为蓝色，0为白色，1为红色
    limits = c(-0.5, 0.7)                         # 设置颜色范围
  ) +
  labs(x = NULL, y = NULL, fill = "Correlation") +
  geom_text(aes(label = rho), color = "black",size=2.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  expand_limits(x = c(0,10)) +  # 扩展 x 轴的范围以显示 p_value
  geom_text(aes(x = 1.5, label = p_value), hjust = 0, nudge_x = 0.3, color = "black")+coord_fixed() # 热图右侧显示 p_value
p4=ggplot(Spatial_TME, aes(x = "", y = reorder(name, rho), fill = rho)) +
  geom_tile(color = "white", width = 1,height = 1) +  # 绘制热图
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),       # 指定颜色
    values = scales::rescale(c(-0.5, 0, 0.7)),    # -1映射为蓝色，0为白色，1为红色
    limits = c(-0.5, 0.7)                         # 设置颜色范围
  ) +
  labs(x = NULL, y = NULL, fill = "Correlation") +
  geom_text(aes(label = rho), color = "black",size=2.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  expand_limits(x = c(0,10)) +  # 扩展 x 轴的范围以显示 p_value
  geom_text(aes(x = 1.5, label = p_value), hjust = 0, nudge_x = 0.3, color = "black")+coord_fixed() # 热图右侧显示 p_value
# 5X5
p1
p2
p3
p4

#FigureS7b
##############ER and Immune#############
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
           "FCGR3A","FCGR3B")
data[,variable]=scale(data[,variable])
data=data[!is.na(data$CAF),]
RNA=data[,c("ER","Arm","Response","TCR_clonality","BCR_clonality",variable)]
RNA <- reshape2::melt(RNA,id.vars=c("ER","Arm","Response"))
#DNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
variable=c("meanHED","TCRA.tcell.fraction.adj")
data[,c("meanHED","TCRA.tcell.fraction.adj")]=scale(data[,c("meanHED","TCRA.tcell.fraction.adj")])
data=data[!is.na(data$TCRA.tcell.fraction.adj),]
DNA=data[,c("ER","Arm","Response","meanHED","TCRA.tcell.fraction.adj")]
DNA <- reshape2::melt(DNA,id.vars=c("ER","Arm","Response"))
#integrative
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
variable=c("Neoantigen")
data[,c("Neoantigen")]=scale(data[,c("Neoantigen")])
data=data[!is.na(data$Neoantigen),]
Inte=data[,c("ER","Arm","Response","Neoantigen")]
Inte<- reshape2::melt(Inte,id.vars=c("ER","Arm","Response"))
#image
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/image_metrics_PREDIX_HER2.rds")
varibale=c("Immune_Cell_prop","Distance_tumor_immune","Cell_Interaction")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin=clin[,c("patientID","ER","Arm","Response")]
image=left_join(data,clin,by="patientID")
image$patientID=NULL
image<- reshape2::melt(image,id.vars=c("ER","Arm","Response"))

d=rbind(RNA,DNA)%>%rbind(Inte)%>%rbind(image)
d$variable=factor(d$variable,levels =c("Neoantigen","meanHED","lohhla","MHC.I_19272155","A01","A03","B07","B44",
                                       "B-cells","DC","Macrophages","T-cells","CD8-T-cells","Cytotoxic-cells","NK-cells","CD45","TILs","Immune_Cell_prop",
                                       "TCR_clonality","BCR_clonality","TCRA.tcell.fraction.adj","FCGR3A","FCGR3B",
                                       "Dysfunction","Exclusion","CAF","Neutrophils","Mast-cells","Treg","TAM_M2","Th2 cells",
                                       "Distance_tumor_immune","Cell_Interaction"))

baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source(paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
Fig <-
  ggplot(d,aes(x=ER,y=value,fill=ER))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=5)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_ER(name="ER status")+
  stat_compare_means(aes(group=ER),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.90)+
  labs(y="",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
Fig

# P 15X11
df=d[d$variable%in%c("TCR_clonality","BCR_clonality","TCRA.tcell.fraction.adj","FCGR3A","FCGR3B",
                     "Dysfunction","Exclusion","CAF","Neutrophils","Mast-cells","Treg","TAM_M2","Th2 cells",
                     "Distance_tumor_immune","Cell_Interaction"),]

source(paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
Fig <-
  ggplot(df,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~ER,scales="free_y",nrow=5)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.90)+
  labs(y="",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
Fig



#FigureS7c
##############Immune and survival#############
library(Blasso);library(tidyverse)
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
Immune=data[,c("patientID",variable)]
colnames(Immune)[1]="ID"
clin=data[,c("patientID","EFS.status",'EFS.time')] 
colnames(clin)=c("ID","status","time")
res<-best_predictor_cox(target_data = clin, 
                        features = Immune, 
                        status = "status",
                        time = "time",
                        nfolds = 5,
                        permutation = 1000)
res
library("survival");library(survminer)
# 1. Determine the optimal cutpoint of variables
# TCR_clonality; CAF
res.cut <- surv_cutpoint(data, time = "EFS.time", event = "EFS.status",variables = c("CAF"))
res.cat <- surv_categorize(res.cut)
head(res.cat)
res.cat$ER=data$ER
res.cat$TUMSIZE=data$TUMSIZE
res.cat$ANYNODES=data$ANYNODES
res.cat$Arm=data$Arm
# 4. Fit survival curves and visualize
fit <- survfit(Surv(EFS.time,EFS.status) ~ CAF, data = res.cat)
ggsurvplot(fit,data = res.cat, palette = c("#BC3C29FF","#6F99ADFF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(res.cat$CAF) # 3X2.5
library(tableone)
res.cat$CAF=factor(res.cat$CAF,levels = c("low","high"))
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(CAF)+Arm+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=res.cat) ##DII_density_with_supp
ShowRegTable(cox.test)

res.cut <- surv_cutpoint(data, time = "EFS.time", event = "EFS.status",variables = c("TCR_clonality"))
res.cat <- surv_categorize(res.cut)
head(res.cat)
res.cat$ER=data$ER
res.cat$TUMSIZE=data$TUMSIZE
res.cat$ANYNODES=data$ANYNODES
res.cat$Arm=data$Arm
# 4. Fit survival curves and visualize
fit <- survfit(Surv(EFS.time,EFS.status) ~ TCR_clonality, data = res.cat)
ggsurvplot(fit,data = res.cat, palette = c("#BC3C29FF","#6F99ADFF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(res.cat$TCR_clonality) # 3X2.5
library(tableone)
res.cat$TCR_clonality=factor(res.cat$TCR_clonality,levels = c("low","high"))
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(TCR_clonality)+Arm+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=res.cat) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))
#DNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
res.cut <- surv_cutpoint(data, time = "EFS.time", event = "EFS.status",variables = c("TCRA.tcell.fraction.adj"))
res.cat <- surv_categorize(res.cut)
head(res.cat)
res.cat$ER=data$ER
res.cat$TUMSIZE=data$TUMSIZE
res.cat$ANYNODES=data$ANYNODES
res.cat$Arm=data$Arm
res.cat$TCRA=res.cat$TCRA.tcell.fraction.adj
# 4. Fit survival curves and visualize
fit <- survfit(Surv(EFS.time,EFS.status) ~ TCRA, data = res.cat)
ggsurvplot(fit,data = res.cat, palette = c("#BC3C29FF","#6F99ADFF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(res.cat$TCRA) # 3.5X3
library(tableone)
res.cat$TCRA=factor(res.cat$TCRA,levels = c("low","high"))
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(TCRA)+Arm+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=res.cat) ##DII_density_with_supp
ShowRegTable(cox.test)



#FigureS7d
library(sunburstR)
library(htmltools)
library(d3r)
library(data.table);library(stringr);library(tidyverse);library(readxl)
Genotype=fread("E:/Projects/PREDIX_HER2/pVACseq/data/Optitype_HLA_PREDIX_HER2.txt")%>%as.data.frame()
Genotype$patientID=substr(Genotype$sampleID,9,12)%>%as.character()
ref=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Resource/Supertype_HAL.xls")
ref$Allele <- sub("(A\\*\\d{2})(\\d{2})", "\\1:\\2", ref$Allele)
ref$Allele <- sub("(B\\*\\d{2})(\\d{2})", "\\1:\\2", ref$Allele)
Genotype$Allele=Genotype$A1
Genotype=left_join(Genotype,ref,by="Allele")
Genotype$SupertypeA1=Genotype$Supertype;Genotype$Supertype=NULL
Genotype$Allele=Genotype$A2
Genotype=left_join(Genotype,ref,by="Allele")
Genotype$SupertypeA2=Genotype$Supertype;Genotype$Supertype=NULL
Genotype$SupertypeA1[Genotype$SupertypeA1%in%c("A01 A24","A01 A03")]="Unclassified"
Genotype$SupertypeA2[Genotype$SupertypeA2%in%c("A01 A24","A01 A03")]="Unclassified"
Genotype$Allele=Genotype$B1
Genotype=left_join(Genotype,ref,by="Allele")
Genotype$SupertypeB1=Genotype$Supertype;Genotype$Supertype=NULL
Genotype$Allele=Genotype$B2
Genotype=left_join(Genotype,ref,by="Allele")
Genotype$SupertypeB2=Genotype$Supertype;Genotype$Supertype=NULL
table(Genotype$SupertypeB1)
table(Genotype$SupertypeB2) # A01 A03
table(Genotype$SupertypeA1)
table(Genotype$SupertypeA2) # A01 A03
tab=table(Genotype[,c("SupertypeA1")],Genotype[,c("SupertypeA2")])%>%as.data.frame()
colnames(tab)=c("Allele1","Allele2","size")
tree <- d3_nest(tab, value_cols = "size")
sunburst(tree, width="100%", height=400)

tab=table(Genotype[,c("SupertypeB1")],Genotype[,c("SupertypeB2")])%>%as.data.frame()
colnames(tab)=c("Allele1","Allele2","size")
tree <- d3_nest(tab, value_cols = "size")
sunburst(tree, width="100%", height=400)

#FigureS7e
library(data.table);library(tidyverse);library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$patientID=as.character(clin$patientID)
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
meta$patientID=as.character(meta$patientID)
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
transport=c("PTPRC","CD3E","GZMA","CD163","CD84","S100A12","CPA3")
MS=MS[transport,]%>%t()%>%as.data.frame()
MS$patientID=row.names(MS)
df=left_join(MS,clin,by="patientID")
df=df[,c(transport,"Response","Arm")]
d <- reshape2::melt(df,id.vars=c("Response","Arm"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
p <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=2)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=5,
                     color="black", label.y.npc = 0.97)+
  labs(y="Normalized protein abundance",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
p


#FigureS7f
library(ggpubr)
library(data.table)
library(ggplot2)
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
Interact_result=separate(Interact_result,biomarker, into = c("group", "per"), sep = "_per_", remove = FALSE)
Interact_result$per=as.numeric(Interact_result$per)
Interact_result$interaction_coefficient=as.numeric(Interact_result$interaction_coefficient)
Interact_result$interaction_coefficient_LCI=as.numeric(Interact_result$interaction_coefficient_LCI)
Interact_result$interaction_coefficient_UCI=as.numeric(Interact_result$interaction_coefficient_UCI)
# ABC_transporter_per_70, Exosome_per_45, Hypoxia_per_80, pik3ca_sig_per_80

pd <- position_dodge(3)
ggplot(Interact_result, 
       aes(x = per, 
           y = interaction_coefficient, group=group, color=group)) +
  geom_point(position=pd, 
             size=3) +
  geom_line(position=pd, 
            size = 1) +
  geom_errorbar(aes(ymin = interaction_coefficient_LCI, 
                    ymax = interaction_coefficient_UCI), 
                width = .1, 
                position=pd, 
                size=1) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal() +
  labs(title = "",
       subtitle = "(mean +/- standard error)",
       x = "", 
       y = "",
       color = "Immunogenomic Biomarker")




















