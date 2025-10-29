######Fig3a
require(openxlsx);library(ggplot2);library(data.table);library(tidyverse);library(glmmSeq)
driverGenes<- scan("E:/Projects/PREDIX_HER2/Multimodal/Resource/breast-cancer-driver-genes.txt", what=character()) #breast cancer
cancergenes<-fread('E:/Projects/PREDIX_HER2/Multimodal/Resource/Census_allMon Nov 28 13_48_06 2022.csv')
cancergenes=cancergenes$`Gene Symbol`
driverGenes=c(driverGenes,cancergenes)
res_TDM1=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=4)
res_TDM1$driver="No"
res_TDM1$driver[res_TDM1$gene%in%driverGenes]="Yes"
res_DHP=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=6)
nrow(res_TDM1[!is.na(res_TDM1$padj)&res_TDM1$padj<0.05&res_TDM1$log2FoldChange>0.5,])
nrow(res_TDM1[!is.na(res_TDM1$padj)&res_TDM1$padj<0.05&res_TDM1$log2FoldChange<(-0.5),])
nrow(res_DHP[!is.na(res_DHP$padj)&res_DHP$padj<0.05&res_DHP$log2FoldChange>0.5,])
nrow(res_DHP[!is.na(res_DHP$padj)&res_DHP$padj<0.05&res_DHP$log2FoldChange<(-0.5),])

genes_to_showname =intersect(union(res_TDM1$gene[abs(res_TDM1$log2FoldChange)>0.5&res_TDM1$padj<0.05],
                                   res_DHP$gene[abs(res_DHP$log2FoldChange)>0.5&res_DHP$padj<0.05]),driverGenes) 
genes_to_showname=union(genes_to_showname,union(res_TDM1$gene[order(res_TDM1$log2FoldChange)][1:5],
                        res_TDM1$gene[order(-res_TDM1$log2FoldChange)][1:5]))
genes_to_showname=union(genes_to_showname,union(res_DHP$gene[order(res_DHP$log2FoldChange)][1:5],
                                                res_DHP$gene[order(-res_DHP$log2FoldChange)][1:5]))
genes_to_showname<- intersect(union(res_TDM1$gene[abs(res_TDM1$log2FoldChange)>0.5&res_TDM1$padj<0.05],
                                                     res_DHP$gene[abs(res_DHP$log2FoldChange)>0.5&res_DHP$padj<0.05]),genes_to_showname) 
genes_to_showname=union(genes_to_showname,c("ABCC12","ABCC3"))
all.equal(res_TDM1$gene,res_DHP$gene)
data=data.frame(ID=res_TDM1$gene,x=res_DHP$log2FoldChange,x_q=res_DHP$padj,y=res_TDM1$log2FoldChange,y_q=res_TDM1$padj)
data[is.na(data)] <- 0.9
genes_to_showname
intersect(genes_to_showname,data$ID)

genes_to_showname=intersect(genes_to_showname,data$ID)
data$combined_sig=0
data$combined_sig[data$x_q<0.05&data$y_q>0.05&abs(data$x)>0.5]=1 # sig only in one arm
data$combined_sig[data$x_q>0.05&data$y_q<0.05&abs(data$y)>0.5]=1 # sig only in one arm
data$combined_sig[data$x_q<0.05&data$y_q<0.05&abs(data$x)>0.5&abs(data$y)>0.5]=2 # sig for both arm
table(data$combined_sig)
data$group[data$combined_sig==0]="Notsig"
data$group[data$x_q<0.05&data$combined_sig==1]="unique DEG in DHP arm"
data$group[data$y_q<0.05&data$combined_sig==1]="unique DEG in T-DM1 arm"
data$group[data$combined_sig==2&data$x*data$y<0]="opposite DEG"
data$group[data$combined_sig==2&data$x*data$y>0]="shared DEG"
table(data$group)
range(data$x);range(data$y)
data$gene_label <- ""
row.names(data)=data$ID
data[genes_to_showname, ]$gene_label <- genes_to_showname
colours = c('grey', 'green3', 'gold3', 'blue')
fig3a=data

annot <- lapply(genes_to_showname, function(i) {
  row <- data[i, ]
  x <- row$x
  y <- row$y
  z <- sqrt(x^2 + y^2)
  list(x = x, y = y,
       text = i, textangle = 0, ax = x/z*75, ay = -y/z*75,
       font = list(color = "black", size =11),
       arrowcolor = "black", arrowwidth = 0.5, arrowhead = 0, arrowsize = 1.5,
       xanchor = "auto", yanchor = "auto")
})
library(plotly)
p <- plot_ly(data = data, x = ~x, y = ~y, type = 'scatter', 
             mode = 'markers',
             color = ~group, colors = colours,
             marker = list(size = 7, 
                           line = list(width = 0.25, color = 'white')),
             text = data$gene_label, hoverinfo = 'text') %>%
  layout(annotations = annot,
         xaxis = list(title = "DHP Arm logFC(pCR vs RD)",
                      color = 'black'),
         yaxis = list(title = "T-DM1 Arm logFC(pCR vs RD)",
                      color = 'black'),
         font = list(size = 12),
         legend = list(x = 0, y = 1, font = list(color = 'black'))) %>%
  config(edits = list(annotationPosition = FALSE,
                      annotationTail = TRUE,
                      annotationText = TRUE),
         toImageButtonOptions = list(format = "svg"))
p

######Fig3b
library(data.table);library(tidyverse);library(readxl)
MS=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_HarmonizR.tsv")
N_NA=colSums(is.na(MS[,2:ncol(MS)]))
quantile(N_NA)
#load mRNA-protein correlation results
x=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/mRNA_protein_cor.xlsx")
res_TDM1=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=4)
gene1=res_TDM1$gene[!is.na(res_TDM1$padj)&res_TDM1$padj<0.05&abs(res_TDM1$log2FoldChange)>0.5]
res_DHP=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=6)
gene2=res_DHP$gene[!is.na(res_DHP$padj)&res_DHP$padj<0.05&abs(res_DHP$log2FoldChange)>0.5]
row.names(x)=x$gene
length(intersect(x$gene,c(gene1,gene2)))
x=x[ intersect(x$gene,c(gene1,gene2)),]
fig3b=x
library(ggplot2)
br <- c(seq(min(x$r,na.rm = TRUE)-0.1,-0.000000001,by=0.05),
        seq(0,max(x$r,na.rm = TRUE),by=0.05))
max_count <- max(table(cut(x$r,breaks = br)))
positive_cor_ratio <- sum(x$r>0,na.rm = TRUE)/nrow(x)

sig_positive_cor_ratio <- sum(x$FDR <= 0.05,na.rm = TRUE)/nrow(x)

mean_cor <- mean(x$r,na.rm = TRUE)
median_cor <- median(x$r,na.rm = TRUE)
n_pairs <- nrow(x)

positive_cor_ratio <- sprintf("%.2f%%",100*positive_cor_ratio)
sig_positive_cor_ratio <- sprintf("%.2f%%",100*sig_positive_cor_ratio) #%>% paste("*\\'%\\'~",sep = "")

x %>% mutate(col=ifelse(r>0,"#BC3C29FF","#2D6DB1")) %>%
  ggplot(aes(x=r,fill=col)) +
  geom_histogram(breaks=br,color="black",size=0.1)+
  xlab("Spearman Correlation Coefficient")+
  ylab("Frequency")+
  ggpubr::theme_pubr()+
  theme(legend.position = "none")+
  scale_fill_manual(breaks = c("#BC3C29FF", "#2D6DB1"),
                    values=c("#BC3C29FF", "#2D6DB1"))+
  geom_vline(xintercept = mean_cor,linetype=2)+
  annotate("text",x=mean_cor,y=1.05*max_count,label=paste("Median = ",sprintf("%.2f",median_cor),sep = ""),
           hjust=-0.05,size=3.5)+
  annotate("text",x=min(br)+0.01,y=1.05*max_count,
           label=paste0(n_pairs," pairs\n",
                        positive_cor_ratio," positive correlation\n",
                        sig_positive_cor_ratio," significant\npositive correlation\n(adjusted P <= 0.05)"),
           vjust=1,hjust=0,size=3.5)


#####Fig3c
library(fgsea)
Hallmark<- gmtPathways("E:/Projects/PREDIX_HER2/Multimodal/Resource/h.all.v2023.2.Hs.symbols.gmt")
geneList=Hallmark
# fgsea - get Normalised Enrichment Score
# For DHP
df=res_DHP[!is.na(res_DHP$padj),]
df$enrichment_score <-df$log2FoldChange
df <- df[order(df$enrichment_score, decreasing = TRUE),]
ranks <- df$enrichment_score
names(ranks) <- df$gene
q <- fgsea(pathways = geneList, stats = ranks,minSize=15, maxSize=500, nperm=100000)
gse_dhp=as.data.frame(q)
gse_dhp$pathway <- gsub("HALLMARK_","",gse_dhp$pathway)
row.names(gse_dhp)=gse_dhp$pathway
# For T-DM1
df=res_TDM1[!is.na(res_TDM1$padj),]
df$enrichment_score <-df$log2FoldChange
df <- df[order(df$enrichment_score, decreasing = TRUE),]
ranks <- df$enrichment_score
names(ranks) <- df$gene
q <- fgsea(pathways = geneList, stats = ranks,minSize=15, maxSize=500, nperm=100000)
gse_tdm1=as.data.frame(q)
gse_tdm1$pathway <- gsub("HALLMARK_","",gse_tdm1$pathway)
row.names(gse_tdm1)=gse_tdm1$pathway
# integrate
term=union(gse_dhp$pathway[gse_dhp$padj<0.05&abs(gse_dhp$NES)>1],gse_tdm1$pathway[gse_tdm1$padj<0.05&abs(gse_tdm1$NES)>1])
gse_dhp=gse_dhp[term,]
gse_tdm1=gse_tdm1[term,]
fig3c=cbind(gse_dhp,gse_tdm1)

# Define the desired order for 'pathway'
gse_dhp=gse_dhp[order(gse_dhp$NES),]
desired_order <-gse_dhp$pathway  # Define your pathway names in the desired order
# Reorder 'pathway' based on the desired order
gse_dhp$pathway <- factor(gse_dhp$pathway, levels = desired_order)
gse_tdm1$pathway <- factor(gse_tdm1$pathway, levels = desired_order)
# ploting
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
dhp <- 
  ggplot(gse_dhp,aes(x=pathway,y=(NES),color=NES))+
  geom_point(aes(size= -log10(padj)))+
  geom_hline(yintercept=0,color="black", linetype="dotted")+
  labs(x="",y="NES")+
  coord_flip()+
  scale_colour_gradient2(low = "#375E97",mid = "white", high = "#FB6542")+
  scale_size_continuous(name=expression(italic(-log[10]~padj)),breaks = c(1:3))+
  theme_manuscript(base_size = figure_font_size)+
  theme(axis.text.y = element_text(face="italic"),
        plot.margin = unit(c(1.1,1.5,0,0.5), "lines"))
tdm1 <- 
  ggplot(gse_tdm1,aes(x=pathway,y=(NES),color=NES))+
  geom_point(aes(size= -log10(padj)))+
  geom_hline(yintercept=0,color="black", linetype="dotted")+
  labs(x="",y="NES")+
  coord_flip()+
  scale_colour_gradient2(low = "#375E97",mid = "white", high = "#FB6542")+
  scale_size_continuous(name=expression(italic(-log[10]~padj)),breaks = c(1:3))+
  theme_manuscript(base_size = figure_font_size)+
  theme(axis.text.y = element_text(face="italic"),
        plot.margin = unit(c(1.1,1.5,0,0.5), "lines"))

ggarrange(dhp,tdm1,nrow = 1,
          font.label = list(size = figure_font_size, family="Helvetica"),
          common.legend =T)

# landscape 11X6


####Fig3.d 
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
rna$patientID=as.double(rna$patientID)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
data=left_join(rna,clin,by="patientID")%>%as.data.frame()
variable=c("Taxane_response","pik3ca_sig","ABC_transporter","Apoptosis","Lysosome","Endocytosis","EMT","Exosome",
           "Oxidative_phosphorylation","Purine_metabolism","Citrate_cycles","Glutathione_metabolism","Fatty_acid_metabolism",
           "Glycolysis","Hypoxia")
norm_variable=c("pik3ca_sig","ABC_transporter","Apoptosis","Lysosome","Endocytosis","EMT","Exosome",
                "Oxidative_phosphorylation","Purine_metabolism","Citrate_cycles","Glutathione_metabolism","Fatty_acid_metabolism",
                "Glycolysis","Hypoxia")
data[,norm_variable]=scale(data[,norm_variable]) 

results=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()
fig3d=results
TDM1=results[,c("biomarker","TDM1_OR","TDM1_lr_p")]
TDM1$group="TDM1"
DHP=results[,c("biomarker","DHP_OR","DHP_lr_p")]
DHP$group="DHP"
colnames(TDM1)=c("Signature","OR","Pvalue","group")
colnames(DHP)=c("Signature","OR","Pvalue","group")
df=rbind(TDM1,DHP)
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
df$logOR=log(as.numeric(df$OR))
df$group=factor(df$group,levels = c("TDM1","DHP"))
unique(df$Signature)
df$Signature=factor(df$Signature,levels =c("Hypoxia","Glycolysis","Fatty_acid_metabolism","Glutathione_metabolism",
                                           "Citrate_cycles","Purine_metabolism","Oxidative_phosphorylation","Apoptosis","EMT","Exosome",
                                           "Endocytosis","Lysosome","ABC_transporter","pik3ca_sig","Taxane_response"))
df=df[order(df$Signature),]
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=df,aes(Signature,logOR,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_modelbarplot(name="Treatment Arm")+
  labs(x="Gene Signature",y="lnOR")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.5,0.5,1), "lines"))+
  coord_flip()+scale_y_continuous(breaks = c(-0.5, 0, 0.5))


#######Fig.3e 
library(data.table)
library(tidyverse)
library(tableone)
library(forestploter)
require(openxlsx)
## continuous variable ##
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
rna=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt')
genomic=inner_join(rna,clin,by="patientID")%>%as.data.frame()
variable=c("Taxane_response","HER2DX_pCR_likelihood_score","pik3ca_sig","ABC_transporter","Exosome","Hypoxia","EMT")
genomic=slice_metric(genomic,variable,20,80,5)
## re-run the batch logistic
selected_columns <- grep("_per_", colnames(genomic), value = TRUE)
Interact_result=Logistic_batch_adjER(genomic,"pCR","Arm",selected_columns ,"ER")%>%as.data.frame()
genomic$Arm=factor(genomic$Arm,levels = c("DHP","T-DM1"))
genomic=as.data.frame(genomic)
res=Logistic_batch_continuous_subgroup(genomic,selected_columns)%>%as.data.frame()

res_continuous=res[res$biomarker%in%c("pik3ca_sig_per_80","ABC_transporter_per_65","Exosome_per_45","Hypoxia_per_80"),]
Interact_continuous=Interact_result[Interact_result$biomarker%in%c("pik3ca_sig_per_80","ABC_transporter_per_65","Exosome_per_45","Hypoxia_per_80"),]
continuous=cbind(res_continuous,Interact_continuous)
# merge
library(openxlsx)
list_of_datasets <- list("rna" = continuous)
write.xlsx(list_of_datasets, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure3/Predictive_biomarker.xlsx")
table(genomic$Hypoxia_per_80,genomic$Arm)
#whole <- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = genomic[genomic$pik3ca_sig_per_80=="Low",])
#ShowRegTable(whole)
#whole <- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = genomic[genomic$pik3ca_sig_per_80=="High",])
#ShowRegTable(whole)
# forest plot #
library(data.table)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure3/Subgroup_Forestplot.csv")
df$DHP <- ifelse(is.na(df$DHP), "", df$DHP)
df$`T-DM1` <- ifelse(is.na(df$`T-DM1`), "", df$`T-DM1`)
df$`P for interaction` <- ifelse(is.na(df$`P for interaction`), "", df$`P for interaction`)
df$` ` <- paste(rep(" ", 20), collapse = " ")
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

#7X5

###Fig3f
###############################################
#Veen correlation between cna-rna cna-protein#
###############################################
library(readxl);library(data.table);library(tidyverse)
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
clinical_anno=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%dplyr::filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample

AMP=c("1p31.3","1q32.2","3p11.1","5p11","6q21","7q11.21","8p11.23","8q23.3","10q22.3","11q13.3",
      "12p11.1","12q15","17q12","19q11","20q13.2")
#DEL=c("1p36.32","1p21.1","1p12","1q21.1","1q21.1","2p11.2","2q13","3p21.2","4p16.1","4q12","4q35.2",  
#      "5q13.2","5q35.2","6q27","7p11.2","7q11.21","7q35","8p23.1","8p21.2","8q24.3","9p24.3","9q11",    
#      "9q13","9q21.11","10q11.22","10q11.22","10q23.2","11p15.5","11p15.4","11q14.3","11q25","12q24.33",
#      "13q11","13q14.13","13q34","14q11.2","14q32.11","15q25.2","16p12.3","16p11.2","16q22.2","17p13.1",
#      "17p11.2","17q21.31","18q23","19p13.3","20q11.1","21p11.2","21q22.12","22p11.2","22q13.33"
#)

gene=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_thresholded.by_genes.txt"))
amp=gene[gene$Cytoband%in%AMP,c("Gene Symbol","Locus ID","Cytoband",sampleid)]
amp_gene=amp$`Gene Symbol`
CNA_mRNA=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx",sheet=1)
CNA_protein=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx",sheet=2)
#Located in amplification focal peaks 
amp_gene=intersect(CNA_mRNA$gene,amp_gene)%>%unique()
CNA_mRNA=CNA_mRNA$gene[CNA_mRNA$gene%in%amp_gene&CNA_mRNA$cis_group=="cis_cor"]
CNA_protein=CNA_protein$gene[CNA_protein$gene%in%amp_gene&CNA_protein$cis_group=="cis_cor"]

x=list(CNA_mRNA_correlation=CNA_mRNA,
       CNA_protein_correlation=CNA_protein) 

library(ggplot2)
library(ggVennDiagram)
ggVennDiagram(x) + scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")

###Fig3g
###############################################
#               CNA drivers (AMP)             #
###############################################
amp=amp[amp$`Gene Symbol`%in%intersect(CNA_mRNA,CNA_protein),]
# select gain or amp
index <- apply(amp[,4:ncol(amp)], 1, function(row) sum(row > 1))>10
amp=amp[index,] 
amp$AMP_rate=apply(amp[,4:ncol(amp)], 1, function(row) sum(row > 1))/179
CNA_mRNA=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx",sheet=1)
row.names(CNA_mRNA)=CNA_mRNA$gene
CNA_mRNA=CNA_mRNA[amp$`Gene Symbol`,]
CNA_protein=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx",sheet=2)
row.names(CNA_protein)=CNA_protein$gene
CNA_protein=CNA_protein[amp$`Gene Symbol`,]
#### add correlation ####
amp$CNA_RNA_cor=CNA_mRNA$r
amp$CNA_protein_cor=CNA_protein$r
#### add lnOR #####
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample
#gene-level
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
CNA=CNA[CNA$`Gene Symbol`%in%amp$`Gene Symbol`,c("Gene Symbol",sampleid)]
row.names(CNA)=CNA$`Gene Symbol`
CNA$`Gene Symbol`=NULL
CNA=as.data.frame(t(CNA))
CNA$patientID=substr(row.names(CNA),1,4)%>%as.integer() 
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
data=left_join(CNA,clin,by='patientID')
# input for 
variable=colnames(CNA)[1:(ncol(CNA)-1)]
##Function without scale
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
CNA_logi=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()
CNA_logi[,2:ncol(CNA_logi)]=apply(CNA_logi[,2:ncol(CNA_logi)],2,as.numeric)
CNA_logi$DHP_OR=CNA_logi$DHP_OR%>%log(base=exp(1))%>%round(digits = 2)
CNA_logi$TDM1_OR=CNA_logi$TDM1_OR%>%log(base=exp(1))%>%round(digits = 2)
#### add logFC RNA and protein #####
library(readxl)
RNA_FC_TDM1=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet = 4)%>%as.data.frame()
RNA_FC_DHP=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet = 6)%>%as.data.frame()
Protein_FC_TDM1=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEqMS_output_TDM1_pCR_RD.tsv")%>%as.data.frame()
Protein_FC_DHP=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEqMS_output_DHP_pCR_RD.tsv")%>%as.data.frame()
row.names(RNA_FC_TDM1)=RNA_FC_TDM1$gene;row.names(RNA_FC_DHP)=RNA_FC_DHP$gene
row.names(Protein_FC_TDM1)=Protein_FC_TDM1$gene;row.names(Protein_FC_DHP)=Protein_FC_DHP$gene
RNA_FC_TDM1=RNA_FC_TDM1[CNA_logi$biomarker,];RNA_FC_DHP=RNA_FC_DHP[CNA_logi$biomarker,];Protein_FC_TDM1=Protein_FC_TDM1[CNA_logi$biomarker,];Protein_FC_DHP=Protein_FC_DHP[CNA_logi$biomarker,]
## visualization ##
library(tidyverse)
library(ggh4x)
library(ggnewscale)
library(ggfun)
df=amp[,c("Gene Symbol","Cytoband","AMP_rate","CNA_RNA_cor","CNA_protein_cor")]
df$AMP_rate=round(100*df$AMP_rate,digits =0)
df$CNA_RNA_cor=round(df$CNA_RNA_cor,digits =2)
df$CNA_protein_cor=round(df$CNA_protein_cor,digits =2)
all.equal(df$`Gene Symbol`,CNA_logi$biomarker,RNA_FC_TDM1$gene,RNA_FC_DHP$gene,Protein_FC_TDM1$gene,Protein_FC_DHP$gene)
df$CNA_DHP_lnOR=CNA_logi$DHP_OR%>%round(digits = 2);df$CNA_DHP_p=CNA_logi$DHP_lr_p%>%as.numeric()%>%round(digits = 2);
df$CNA_TDM1_lnOR=CNA_logi$TDM1_OR%>%round(digits = 2);df$CNA_TDM1_p=CNA_logi$TDM1_lr_p%>%as.numeric()%>%round(digits = 2);df$CNA_interact_p=CNA_logi$interaction_lr_p%>%as.numeric()%>%round(digits = 2);
df$RNA_DHP_logFC=RNA_FC_DHP$log2FoldChange%>%round(digits = 2);df$RNA_DHP_padj=RNA_FC_DHP$padj%>%round(digits = 2)
df$RNA_TDM1_logFC=RNA_FC_TDM1$log2FoldChange%>%round(digits = 2);df$RNA_TDM1_padj=RNA_FC_TDM1$padj%>%round(digits = 2)
df$Prot_DHP_logFC=Protein_FC_DHP$logFC%>%round(digits = 2);df$Prot_DHP_p=Protein_FC_DHP$P.Value%>%round(digits = 2)
df$Prot_TDM1_logFC=Protein_FC_TDM1$logFC%>%round(digits = 2);df$Prot_TDM1_p=Protein_FC_TDM1$P.Value%>%round(digits = 2)
fig3g=df
#### transfer data ####
df2 <- df %>% tidyr::pivot_longer(cols = -c("Gene Symbol","Cytoband","CNA_DHP_p","CNA_TDM1_p","CNA_interact_p",
                                            "RNA_DHP_padj","RNA_TDM1_padj","Prot_DHP_p","Prot_TDM1_p"), 
                                  names_to = "Group", values_to = "Value")
df2$`Gene Symbol`=factor(df2$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
library(ggplot2)
library(cowplot)
df_amp_rate <- subset(df2, Group == "AMP_rate")
df_amp_rate$`Gene Symbol`=factor(df_amp_rate$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
df_cna_rna_cor <- subset(df2, Group == "CNA_RNA_cor")
df_cna_rna_cor$`Gene Symbol`=factor(df_cna_rna_cor$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
df_cna_protein_cor <- subset(df2, Group == "CNA_protein_cor")
df_cna_protein_cor$`Gene Symbol`=factor(df_cna_protein_cor$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
df_FC=subset(df2, Group %in% c("CNA_DHP_lnOR","CNA_TDM1_lnOR","RNA_DHP_logFC","RNA_TDM1_logFC","Prot_DHP_logFC","Prot_TDM1_logFC"))
df_FC$Group=factor(df_FC$Group,levels = c("Prot_TDM1_logFC","Prot_DHP_logFC","RNA_TDM1_logFC","RNA_DHP_logFC","CNA_TDM1_lnOR","CNA_DHP_lnOR"))
df_FC$`Gene Symbol`=factor(df_FC$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
df_FC$Label=NA
df_FC$Label[df_FC$Group=="Prot_TDM1_logFC"] <- with(df_FC[df_FC$Group=="Prot_TDM1_logFC",], ifelse(Prot_TDM1_p < 0.01, paste(Value, "\n***", sep=""),
                                                                                                   ifelse(Prot_TDM1_p < 0.05, paste(Value, "\n**", sep=""),
                                                                                                          ifelse(Prot_TDM1_p < 0.1, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="Prot_DHP_logFC"] <- with(df_FC[df_FC$Group=="Prot_DHP_logFC",], ifelse(Prot_DHP_p < 0.01, paste(Value, "\n***", sep=""),
                                                                                                 ifelse(Prot_DHP_p < 0.05, paste(Value, "\n**", sep=""),
                                                                                                        ifelse(Prot_DHP_p < 0.1, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="RNA_TDM1_logFC"] <- with(df_FC[df_FC$Group=="RNA_TDM1_logFC",], ifelse(RNA_TDM1_padj < 0.01, paste(Value, "\n***", sep=""),
                                                                                                 ifelse(RNA_TDM1_padj < 0.05, paste(Value, "\n**", sep=""),
                                                                                                        ifelse(RNA_TDM1_padj < 0.1, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="RNA_DHP_logFC"] <- with(df_FC[df_FC$Group=="RNA_DHP_logFC",], ifelse(RNA_DHP_padj < 0.01, paste(Value, "\n***", sep=""),
                                                                                               ifelse(RNA_DHP_padj < 0.05, paste(Value, "\n**", sep=""),
                                                                                                      ifelse(RNA_DHP_padj < 0.1, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="CNA_TDM1_lnOR"] <- with(df_FC[df_FC$Group=="CNA_TDM1_lnOR",], ifelse(CNA_TDM1_p < 0.01, paste(Value, "\n***", sep=""),
                                                                                               ifelse(CNA_TDM1_p < 0.05, paste(Value, "\n**", sep=""),
                                                                                                      ifelse(CNA_TDM1_p <= 0.11, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="CNA_DHP_lnOR"] <- with(df_FC[df_FC$Group=="CNA_DHP_lnOR",], ifelse(CNA_DHP_p < 0.01, paste(Value, "\n***", sep=""),
                                                                                             ifelse(CNA_DHP_p < 0.05, paste(Value, "\n**", sep=""),
                                                                                                    ifelse(CNA_DHP_p <= 0.11, paste(Value, "\n*", sep=""), as.character(Value)))))
p1 <- ggplot(df_amp_rate, aes(`Gene Symbol`, Group, fill=Value)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=paste0(Value,"%")), color="black", size=4) +
  scale_fill_gradientn(colours=colorRampPalette(c('#e9e9e9', '#d4d4d4', '#707070'))(111)) +
  labs(x=NULL, y=NULL, fill="AMP rate (%)") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

p2 <- ggplot(df_cna_rna_cor, aes(`Gene Symbol`, Group, fill=Value)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=Value), color="gray2", size=4) +
  scale_fill_gradientn(colours=colorRampPalette(c('#d4c0d9', '#a982b4', '#643d6e'))(111)) +
  labs(x=NULL, y=NULL, fill="CNA-RNA correlation") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

p3 <- ggplot(df_cna_protein_cor, aes(`Gene Symbol`, Group, fill=Value)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=Value), color="gray2", size=4) +
  scale_fill_gradientn(colours=colorRampPalette(c('#FFCDB2', '#D06814'))(111)) +
  labs(x=NULL, y=NULL, fill="CNA-protein correlation") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
library(scales)
p4 <- ggplot(df_FC, aes(`Gene Symbol`, Group, fill=Value)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=Label), color="gray2", size=4) +
  scale_fill_gradientn(
    colours=c('#2D6DB1', 'white', '#DC1623'),
    values=rescale(c(-0.5,  0, 2))
  ) +
  labs(x=NULL, y=NULL, fill="lnOR/logFC") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9, size=9))

legend1 <- get_legend(p1) 
legend2 <- get_legend(p2)
legend3 <- get_legend(p3)
legend4 <- get_legend(p4)
combined_legend <- plot_grid(legend1, legend2, legend3, legend4, ncol = 4, align = 'h')

p1 <- p1 + theme(legend.position="none")
p2 <- p2 + theme(legend.position="none")
p3 <- p3 + theme(legend.position="none")
p4 <- p4 + theme(legend.position="none")
combined_plot <- plot_grid(combined_legend, 
                           plot_grid(p1, p2, p3, p4, ncol = 1, align = "v", rel_heights = c(1, 1, 1, 4)),
                           ncol = 1, rel_heights = c(0.2, 0.8))


combined_plot

# 12X7.5

###Fig3h
# early endosome：RAB5C, EEA1
# late endosome: RAB7A, RILP 
# Lysosomal  SLC12A2 CTSL  SLC46A3

# Recycling endosome: ARL1  RAB11B
# ABC transporter ABCC12
# Exosome VAMP3 FLOT1 
library(readr);library(tidyverse);library(ggplot2);library(data.table);library(circlize)
gene=c("ARL1","RAB11B","VAMP3","FLOT1","ABCC12","RAB5C","EEA1","RAB7A","RILP","SLC12A2","CTSL","ABCC1")
res_TDM1=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet="DEG_TDM1")%>%as.data.frame()
row.names(res_TDM1)=res_TDM1$gene;res_TDM1=res_TDM1[gene,c("gene","log2FoldChange","pvalue")];colnames(res_TDM1)=paste0(colnames(res_TDM1),"_rna_tdm1") 
res_DHP=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet="DEG_DHP")%>%as.data.frame()
row.names(res_DHP)=res_DHP$gene;res_DHP=res_DHP[gene,c("gene","log2FoldChange","pvalue")];colnames(res_DHP)=paste0(colnames(res_DHP),"_rna_dhp") 
prot_TDM1=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEqMS_output_TDM1_pCR_RD.tsv")%>%as.data.frame()
row.names(prot_TDM1)=prot_TDM1$gene;prot_TDM1=prot_TDM1[gene,c("gene","logFC","P.Value")];colnames(prot_TDM1)=paste0(colnames(prot_TDM1),"_prot_tdm1") 
prot_DHP=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEqMS_output_DHP_pCR_RD.tsv")%>%as.data.frame()
row.names(prot_DHP)=prot_DHP$gene;prot_DHP=prot_DHP[gene,c("gene","logFC","P.Value")];colnames(prot_DHP)=paste0(colnames(prot_DHP),"_prot_dhp") 
#col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c("#00FF00", "#008000", "#000000", "#800000", "#FF0000"))
col_fun = colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#4575b4", "#91bfdb", "#f7f7f7", "#fc8d59", "#d73027")
)
data=cbind(res_TDM1,res_DHP,prot_TDM1,prot_DHP)
fig3h=data
# logFC 
mat = data[, c("log2FoldChange_rna_dhp","logFC_prot_dhp","log2FoldChange_rna_tdm1","logFC_prot_tdm1")]
rownames(mat) = data$gene_rna_tdm1   
mat = as.matrix(mat)
colnames(mat) = c("RNA_DHP","Prot_DHP", "RNA_TDM1", "Prot_TDM1")

library(ComplexHeatmap)
ht = Heatmap(
  mat,
  name = "logFC",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "logFC", at = c(-1, -0.5, 0, 0.5, 1)),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), 
              x, y, gp = gpar(fontsize = 12, col = "black"))
  }
)


ht

# Portrait 4.5 X 5.5 25%

library(openxlsx)
data_list=list("fig3a"=fig3a%>%as.data.frame(),
               "fig3b"=fig3b%>%as.data.frame(),
               "fig3c"=fig3c%>%as.data.frame(),
               "fig3d"=fig3d%>%as.data.frame(),
               "fig3g"=fig3g%>%as.data.frame(),
               "fig3h"=fig3h%>%as.data.frame())
openxlsx::write.xlsx(data_list,file='E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Figure3.xlsx')








