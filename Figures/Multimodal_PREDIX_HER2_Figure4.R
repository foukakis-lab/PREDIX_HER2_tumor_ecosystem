# Fig4a
library(networkD3);library(tidyr);library(tibble);library(dplyr)
library(UpSetR);library(ggplot2);library(plyr);library(gridExtra)
library(ggpubr);library(data.table);library(tableone)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(df,clin,by="patientID")%>%as.data.frame()
table(df$sspbc.subtype,df$Response)
fig4a=df[,c("patientID","Arm","sspbc.subtype","Response")]

df$sspbc.bin="Non-HER2"
df$sspbc.bin[df$sspbc.subtype=="Her2"]="HER2-enriched"
df$sspbc.bin=factor(df$sspbc.bin,levels = c("Non-HER2","HER2-enriched"))
res<- glm(as.numeric(pCR) ~ sspbc.bin+ER, family = "binomial", data = df[df$Arm=="T-DM1",])
ShowRegTable(res)
table(df[df$Arm=="T-DM1","Response"],df[df$Arm=="T-DM1","sspbc.bin"])
table(df[df$Arm=="DHP","Response"],df[df$Arm=="DHP","sspbc.bin"])

res<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = df[df$sspbc.subtype=="Her2",])
ShowRegTable(res)
table(df$sspbc.subtype)

table(df$sspbc.subtype,df$Response)
df1=df[,c("Arm","sspbc.subtype")];colnames(df1)=c('source','target')
df2=df[,c("sspbc.subtype","Response")];colnames(df2)=c('source','target')
links=rbind(df1,df2)
links$color=links$source

my_color <- 'd3.scaleOrdinal() .domain(["DHP","T-DM1","LumA","LumB","Her2","Basal","RD","pCR"]) .range(["#8491B4FF","#91D1C2FF","#1f78b4","#a6cee3","#fb9a99","#e31a1c","#fdb462","#00A087FF"])'

links$source <- as.character(links$source)
links$target<- as.character(links$target)
table(links$target)
nodes <- data.frame(name = unique(c(links$source,links$target)))
nodes$color=nodes$name

links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1
links$value <- 1 # add also a value

sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',colourScale=my_color,LinkGroup = NULL,
              Target = 'target', Value = 'value', NodeID = 'name',NodeGroup="color",
              height ="600",width ="800",nodeWidth ="50")


# Fig4b HER2-prot-low
library(data.table);library(tidyverse);library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
MS=MS[c("ERBB2","GRB7","MIEN1"),]%>%t()%>%as.data.frame()
MS$patientID=row.names(MS)
MS$HER2_amplicon=rowSums(MS[,1:3])
pid=intersect(MS$patientID,clin$patientID[clin$ISH_HER2_copy<6&clin$HER2neu1=="2+"])
mean=mean(MS[pid,"HER2_amplicon"])
sd=sd(MS[pid,"HER2_amplicon"])
MS$HER2_low="No"
MS$HER2_low[MS$patientID%in%clin$patientID[clin$ISH_HER2_copy<6&clin$HER2neu1=="2+"]]="Yes"
table(MS$HER2_low)
MS$Zscore=(MS$HER2_amplicon-mean)/sd
MS$HER2_prot=NA
MS$HER2_prot[MS$HER2_low=="Yes"]="HER2_low"
MS$HER2_prot[is.na(MS$HER2_prot)&MS$Zscore<2]="Pseudo_HER2"
MS$HER2_prot[is.na(MS$HER2_prot)]="HER2_pos"
table(MS$HER2_prot)
meta=meta[,c("sampleID","patientID")]
MS$patientID=as.integer(MS$patientID)
MS=left_join(meta,MS,by="patientID")
MS$sampleID=NULL
fig4b=MS
#saveRDS(MS,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_MS.rds")
gghistogram(MS, x = "Zscore", bins = 41, rug = TRUE,
            color = "HER2_low", fill = "HER2_low",
            palette = c("darkgreen","darkorchid4")) +
  geom_vline(xintercept = c(-2,0,2,4,6,8), linetype = "dashed", color = "black", size = 0.2)+
  scale_x_continuous(breaks = seq(-2, 10, by = 2)) +
  scale_y_continuous(breaks = seq(2, 10, by = 2)) +
  # Add density plot for HER2_low == "Yes"
  geom_density(data = subset(MS, HER2_low == "Yes"), aes(x = Zscore, y = ..density.. *15), 
               fill = "darkorchid4", alpha = 0.3, size = 0.5) 
# 5X3


blue2red <- colorRampPalette(c("midnightblue","#003BC0","#6689d9","white","#eb8694","#D80D29","firebrick4"))
green2purple <- colorRampPalette(c("darkgreen","white","darkorchid4"))

grey2purple <- colorRampPalette(c("grey80","#8ecee2","#A074B6","darkorchid4"))
####Prot and pCR ####
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
data=left_join(clin,MS,by="patientID")
table(data$Response,data$HER2_prot,data$Arm)
chisq.test(data$Response[data$HER2neu1=="2+"],data$Arm[data$HER2neu1=="2+"])
chisq.test(data$Response[data$HER2neu1=="3+"],data$Arm[data$HER2neu1=="3+"])


# Fig4c Heatmap HER2 class
library(tableone);library(data.table);library(tidyverse)
# protein heatmap
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
protein=c("ERBB2","GRB7","MIEN1","AKT1","AKT2","MAP2K1","MAPK1","PTEN","TOP2A")
protein=MS[protein,]%>%t()%>%as.data.frame()
protein$patientID=row.names(protein)%>%as.double()
# HER2 protein-low
MS=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_MS.rds")
MS$HER2_prot[MS$HER2_prot=="Pseudo_HER2"]="HER2_low"
MS[,c("ERBB2","GRB7","MIEN1")]=NULL
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$HER2neu1[clin$ISH_HER2_copy<6]
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
rna$patientID=as.double(rna$patientID)
mut=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
mut=mut[,c("patientID","coding_mutation_ERBB2")]
mut$patientID=as.double(mut$patientID)
mut$ERBB2_mut[mut$coding_mutation_ERBB2==1]="Mut"
mut$ERBB2_mut[mut$coding_mutation_ERBB2==0]="Wild-type"
mut$coding_mutation_ERBB2=NULL
# CNA CUTseq
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clin$patientID,]
sampleid=seg_count$sample
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
amp=c("ERBB2")
CNA=CNA[CNA$`Gene Symbol`%in%amp,c("Gene Symbol",sampleid)]
row.names(CNA)=CNA$`Gene Symbol`
CNA$`Gene Symbol`=NULL
CNA=as.data.frame(t(CNA))
CNA$patientID=substr(row.names(CNA),1,4)%>%as.integer() 
colnames(CNA)[1]="ERBB2_CN"
# CNA WES
WES_CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_data_by_genes.txt")%>%as.data.frame()
amp=c("ERBB2")
WES_CNA=WES_CNA[WES_CNA$`Gene Symbol`%in%amp,]
row.names(WES_CNA)=WES_CNA$`Gene Symbol`
WES_CNA$`Gene Symbol`=NULL;WES_CNA$`Gene ID`=NULL;WES_CNA$Cytoband=NULL
WES_CNA=as.data.frame(t(WES_CNA))
WES_CNA$patientID=substr(row.names(WES_CNA),9,12)%>%as.integer() 
colnames(WES_CNA)[1]="ERBB2_CN"
pid=setdiff(WES_CNA$patientID,CNA$patientID)
WES_CNA=WES_CNA[WES_CNA$patientID%in%pid,]
# merge CNA
CNA=rbind(CNA,WES_CNA)
data=left_join(MS,clin,by="patientID")%>%left_join(rna,by="patientID")%>%left_join(CNA,by="patientID")%>%
  left_join(mut,by="patientID")%>%left_join(protein,by="patientID")
data$ERBB2_amp=NA
data$ERBB2_amp[data$ERBB2_CN<2]="No"
data$ERBB2_amp[data$ERBB2_CN>2]="Yes"
table(data$ERBB2_amp)
table(data$ERBB2_amp,data$HER2_prot)
table(data$Response,data$HER2_prot)
table(data$sspbc.subtype,data$HER2_prot)
data$ERBB2_PG="Negative"
data$ERBB2_PG[data$HER2_prot=='HER2_pos'&data$ERBB2_amp=="Yes"]="Positive"
table(data$ERBB2_PG)
data$ERBB2_class=NA
data$ERBB2_class[(data$HER2_prot!='HER2_pos'|data$ERBB2_amp!="Yes")&data$sspbc.subtype=="Her2"]="Her2E ERBB2 PG-"
data$ERBB2_class[(data$HER2_prot!='HER2_pos'|data$ERBB2_amp!="Yes")&data$sspbc.subtype!="Her2"]="Other ERBB2 PG-"
data$ERBB2_class[data$HER2_prot=='HER2_pos'&data$ERBB2_amp=="Yes"&data$sspbc.subtype=="Her2"]="Her2E ERBB2 PG+"
data$ERBB2_class[data$HER2_prot=='HER2_pos'&data$ERBB2_amp=="Yes"&data$sspbc.subtype!="Her2"]="Other ERBB2 PG+"
data$ERBB2_class=factor(data$ERBB2_class,levels = c("Her2E ERBB2 PG-","Other ERBB2 PG-","Her2E ERBB2 PG+","Other ERBB2 PG+"))
table(data$ERBB2_class,data$Response)
table(is.na(data$ERBB2_class))
#saveRDS(data,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
library(data.table);library(tidyverse);library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
library(ComplexHeatmap);library(circlize);library(RColorBrewer);library(ggstatsplot)
df=data[!is.na(data$ERBB2_class),]
df=df[order(df$ERBB2_class),]
row.names(df)=df$patientID
protein=c("ERBB2","GRB7","MIEN1","AKT1","AKT2","MAP2K1","MAPK1","PTEN","TOP2A")
df=df[,c("ERBB2_class","Arm","Response","EFS.status","sspbc.subtype","HER2_prot","ERBB2_PG","ER","HER2neu1","ERBB2_amp","ERBB2_mut",protein)]
#### HER2DX group
s1 = as.matrix(t(df[df$ERBB2_class=="Her2E ERBB2 PG-",protein]))   # subtype1
S1_meta=df[colnames(s1),]
s2 = as.matrix(t(df[df$ERBB2_class=="Other ERBB2 PG-",protein]))   # subtype1
S2_meta=df[colnames(s2),]
s3 = as.matrix(t(df[df$ERBB2_class=="Her2E ERBB2 PG+",protein]))   # subtype1
S3_meta=df[colnames(s3),]
s4 = as.matrix(t(df[df$ERBB2_class=="Other ERBB2 PG+",protein]))   # subtype1
S4_meta=df[colnames(s4),]

pheno=S1_meta
ha1=HeatmapAnnotation(ERBB2_class=pheno$ERBB2_class,
                      Intrinsic_subtype=pheno$sspbc.subtype,
                      ERBB2_PG=pheno$ERBB2_PG,
                      HER2_prot=pheno$HER2_prot,
                      ERBB2_amp=pheno$ERBB2_amp,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      col=list(ERBB2_class=c("Her2E ERBB2 PG-"="#F6BD60","Other ERBB2 PG-"="#F7EDE2","Her2E ERBB2 PG+"="#F5CAC3","Other ERBB2 PG+"="#84A59D"),
                               Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               Intrinsic_subtype= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               HER2_prot=c("HER2_low"="#E8E8E8","HER2_pos"="black"),
                               ERBB2_PG=c("Negative"="#E8E8E8","Positive"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ERBB2_amp=c("No"="#E8E8E8","Yes"="black")),
                      height = unit(0.5, "cm"),na_col="#808080",show_annotation_name = FALSE,
                      border=F)

pheno=S2_meta
ha2=HeatmapAnnotation(ERBB2_class=pheno$ERBB2_class,
                      Intrinsic_subtype=pheno$sspbc.subtype,
                      ERBB2_PG=pheno$ERBB2_PG,
                      HER2_prot=pheno$HER2_prot,
                      ERBB2_amp=pheno$ERBB2_amp,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      col=list(ERBB2_class=c("Her2E ERBB2 PG-"="#F6BD60","Other ERBB2 PG-"="#F7EDE2","Her2E ERBB2 PG+"="#F5CAC3","Other ERBB2 PG+"="#84A59D"),
                               Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               Intrinsic_subtype= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               HER2_prot=c("HER2_low"="#E8E8E8","HER2_pos"="black"),
                               ERBB2_PG=c("Negative"="#E8E8E8","Positive"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ERBB2_amp=c("No"="#E8E8E8","Yes"="black")),
                      height = unit(0.5, "cm"),na_col="#808080",show_annotation_name = FALSE,
                      border=F)

pheno=S3_meta
ha3=HeatmapAnnotation(ERBB2_class=pheno$ERBB2_class,
                      Intrinsic_subtype=pheno$sspbc.subtype,
                      ERBB2_PG=pheno$ERBB2_PG,
                      HER2_prot=pheno$HER2_prot,
                      ERBB2_amp=pheno$ERBB2_amp,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      col=list(ERBB2_class=c("Her2E ERBB2 PG-"="#F6BD60","Other ERBB2 PG-"="#F7EDE2","Her2E ERBB2 PG+"="#F5CAC3","Other ERBB2 PG+"="#84A59D"),
                               Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               Intrinsic_subtype= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               HER2_prot=c("HER2_low"="#E8E8E8","HER2_pos"="black"),
                               ERBB2_PG=c("Negative"="#E8E8E8","Positive"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ERBB2_amp=c("No"="#E8E8E8","Yes"="black")),
                      height = unit(0.5, "cm"),na_col="#808080",show_annotation_name = FALSE,
                      border=F)

pheno=S4_meta
ha4=HeatmapAnnotation(ERBB2_class=pheno$ERBB2_class,
                      Intrinsic_subtype=pheno$sspbc.subtype,
                      ERBB2_PG=pheno$ERBB2_PG,
                      HER2_prot=pheno$HER2_prot,
                      ERBB2_amp=pheno$ERBB2_amp,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      col=list(ERBB2_class=c("Her2E ERBB2 PG-"="#F6BD60","Other ERBB2 PG-"="#F7EDE2","Her2E ERBB2 PG+"="#F5CAC3","Other ERBB2 PG+"="#84A59D"),
                               Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               Intrinsic_subtype= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               HER2_prot=c("HER2_low"="#E8E8E8","HER2_pos"="black"),
                               ERBB2_PG=c("Negative"="#E8E8E8","Positive"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ERBB2_amp=c("No"="#E8E8E8","Yes"="black")),
                      height = unit(0.5, "cm"),na_col="#808080",show_annotation_name = T,
                      border=F)
#### Plot heatmap
col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), c("midnightblue","#003BC0","#6689d9","white","#eb8694","#D80D29","firebrick4"))

ht_list = Heatmap(as.matrix(s1), col = col_fun, name = "Protein abundance",
                  cluster_rows =T,show_row_names = FALSE,cluster_columns = F,
                  show_row_dend = T, show_column_dend = F,
                  show_column_names = F,
                  top_annotation = ha1,
                  row_title_gp = gpar(col = "#FFFFFF00"), width = unit(2.1, "cm"),height= unit(5, "cm"))+
  Heatmap(as.matrix(s2),col=col_fun, show_column_names = F, cluster_columns = F,
          show_column_dend = F,top_annotation=ha2,
          show_heatmap_legend = FALSE, width = unit(3.2, "cm"),height= unit(5, "cm"))+
  Heatmap(as.matrix(s3),col=col_fun, show_column_names = F, cluster_columns = F,
          show_column_dend = F,top_annotation=ha3,row_names_gp = gpar(fontsize = 2),
          show_heatmap_legend = FALSE, width = unit(5.9, "cm"),height= unit(5, "cm"))+
  Heatmap(as.matrix(s4),col=col_fun, show_column_names = F, cluster_columns = F,
          show_column_dend = F,top_annotation=ha4,row_names_gp = gpar(fontsize = 8),
          show_heatmap_legend = FALSE, width = unit(2, "cm"),height= unit(5, "cm"))

p0=draw(ht_list, row_title = paste0("protein abundance"),
        annotation_legend_side = "left", heatmap_legend_side = "left") 
p0

# Fig4d
library(ggplot2);library(ggpubr);library(dplyr);library(tidyverse)
mydata=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
table(mydata$Response[mydata$Arm=="DHP"],mydata$ERBB2_PG[mydata$Arm=="DHP"])
chisq.test(mydata$Response[mydata$Arm=="DHP"],mydata$ERBB2_PG[mydata$Arm=="DHP"]) #0.006
table(mydata$Response[mydata$Arm=="T-DM1"],mydata$ERBB2_PG[mydata$Arm=="T-DM1"])
chisq.test(mydata$Response[mydata$Arm=="T-DM1"],mydata$ERBB2_PG[mydata$Arm=="T-DM1"]) #0.006
mydata$ERBB2_PG=factor(mydata$ERBB2_PG,levels = c("Negative","Positive"))
d=mydata[mydata$Arm=="DHP",c("ERBB2_PG","Response")]
mytable <- d %>%
  count(ERBB2_PG, Response) %>%
  group_by(ERBB2_PG) %>%
  mutate(freq = n / sum(n))
dhp=ggbarplot(mytable, "ERBB2_PG", "freq",
              fill = "Response", color = "Response", 
              palette = c("#fdb462","#00A087FF"))

d=mydata[mydata$Arm=="T-DM1",]
mytable=d %>%
  group_by(ERBB2_PG,Response) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
tdm1=ggbarplot(mytable, "ERBB2_PG", "freq",
               fill = "Response", color = "Response", 
               palette = c("#fdb462","#00A087FF"))


dhp+tdm1
library(tableone)
res<- glm(as.numeric(pCR) ~ ERBB2_PG+ER, family = "binomial", data = mydata[mydata$Arm=="DHP",])
ShowRegTable(res)
res<- glm(as.numeric(pCR) ~ ERBB2_PG+ER, family = "binomial", data = mydata[mydata$Arm=="T-DM1",])
ShowRegTable(res)



# Fig4e
library(tidyverse);library(tableone);library(data.table);library(lmtest);library(forestploter)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
data$HER2_prot[data$HER2_prot!="HER2_pos"]="Low"
data$HER2_prot[data$HER2_prot=="HER2_pos"]="High"

table(data$HER2_prot)
table(data$HER2_prot,data$Response,data$Arm)
whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$HER2_prot=="Low",])
ShowRegTable(whole)
whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$HER2_prot=="High",])
ShowRegTable(whole)
interaction_1<- glm(as.numeric(pCR) ~ HER2_prot+ER+Arm, family = "binomial", data = data)
interaction_2<- glm(as.numeric(pCR) ~ HER2_prot+ER+Arm+Arm*HER2_prot, family = "binomial", data = data)
lrtest(interaction_1,interaction_2)


data$ERBB2_class=factor(data$ERBB2_class,levels = c("Other ERBB2 PG-","Her2E ERBB2 PG-","Her2E ERBB2 PG+","Other ERBB2 PG+"))
table(data$ERBB2_PG)
table(data$ERBB2_PG,data$Arm)
whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$ERBB2_PG=="Negative",])
ShowRegTable(whole)
whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$ERBB2_PG=="Positive",])
ShowRegTable(whole)
interaction_1<- glm(as.numeric(pCR) ~ ERBB2_PG+ER+Arm, family = "binomial", data = data)
interaction_2<- glm(as.numeric(pCR) ~ ERBB2_PG+ER+Arm+Arm*ERBB2_PG, family = "binomial", data = data)
lrtest(interaction_1,interaction_2)

data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
data$ERBB2_class=as.character(data$ERBB2_class)
data$ERBB2_class[data$ERBB2_class%in%c("Other ERBB2 PG-","Her2E ERBB2 PG-","Other ERBB2 PG+")]="Others"

table(data$ERBB2_class)
whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$ERBB2_class=="Her2E ERBB2 PG+",])
ShowRegTable(whole)
table(data[data$ERBB2_class%in%c("Her2E ERBB2 PG+"),"Response"],data[data$ERBB2_class%in%c("Her2E ERBB2 PG+"),"Arm"])

whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$ERBB2_class=="Others",])
ShowRegTable(whole)
table(data[data$ERBB2_class%in%c("Others"),"Response"],data[data$ERBB2_class%in%c("Others"),"Arm"])
interaction_1<- glm(as.numeric(pCR) ~ ERBB2_class+ER+Arm, family = "binomial", data = data)
interaction_2<- glm(as.numeric(pCR) ~ ERBB2_class+ER+Arm+Arm*ERBB2_class, family = "binomial", data = data)
lrtest(interaction_1,interaction_2)

library(tidyverse);library(tableone);library(data.table);library(lmtest);library(forestploter)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure4/Subgroup_Forestplot.csv")
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

# 7.5X5
