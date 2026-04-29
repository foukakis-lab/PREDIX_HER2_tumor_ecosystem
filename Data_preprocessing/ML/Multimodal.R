# prepare data for machine learning
library(data.table)
setwd('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics')
library(data.table);library(tidyverse);library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin=clin[,c('patientID',"Arm","pCR","prolifvalu","TUMSIZE","ANYNODES","ER")] 
clin$prolifvalu[clin$prolifvalu=="NA"]="Unknown"
clin$TUMSIZE[clin$TUMSIZE=="NA"]="Unknown"
clin$ANYNODES[clin$ANYNODES=="NA"]="Unknown"
colnames(clin)[c(2, 4:7)]=paste0("Clin_",colnames(clin)[c(2, 4:7)])

DNA=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
#DNA=na.omit(DNA)
colnames(DNA)
DNA_name=c("ERBB2_CNA","BRCA2_CNA","PIK3CA_CNA","RAB11FIP1_CNA","PPFIA1_CNA","CTTN_CNA",             
       "RPL19_CNA","PPP1R1B_CNA","MIEN1_CNA","GRB7_CNA","coding_mutation_TP53_oncokb",     
       "coding_mutation_PIK3CA_oncokb","coding_mutation_ERBB2_oncokb","coding_mutation_GATA3_oncokb",
       "LOH_Del_burden","CNV_burden","TMB_uniform","TMB_clone","HRD")
DNA <- DNA[, c("patientID",DNA_name)]
mut=c("coding_mutation_TP53_oncokb","coding_mutation_PIK3CA_oncokb","coding_mutation_GATA3_oncokb","coding_mutation_ERBB2_oncokb")
DNA[,mut]=lapply(DNA[, mut], function(x) ifelse(x == 1, "TRUE", "FALSE"))
DNA$TMB_clone[DNA$TMB_clone]=0
library(deconstructSigs);library(BSgenome.Hsapiens.UCSC.hg38);library(data.table)
# https://www.nature.com/articles/s41586-019-1056-z#Sec2
maf=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_vaf0.5_curated.rds")
maf$vaf=maf$t_alt_count/maf$t_depth

#tumor_barcode=freq$Tumor_Sample_Barcode[freq$n>10]
#maf=maf%>%filter(Tumor_Sample_Barcode%in%tumor_barcode)
# Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = maf, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
cosine=t(sigs.input)%>%as.data.frame()
cosine$'Mutation Types'=row.names(cosine)
cosine <- cosine[, c(ncol(cosine), 1:(ncol(cosine)-1))]
write.table(cosine,file='E:/Projects/PREDIX_HER2/Multimodal/Analyses/SigProfilerExtractor/PREDIX_HER2_cosine96.txt',quote = F,row.names =F,sep="\t")

w=lapply(row.names(sigs.input), function(i){
  ## signatures.cosmic signatures.nature2013
  sample_1 = whichSignatures(tumor.ref = sigs.input[,], 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'exome')
  print(i)
  return(sample_1$weights)
})
w=do.call(rbind,w)
df <- w[, colMeans(w == 0) < 0.90]
ncol(df)
df[df==0]=0.001
# normalize columns using the first column (log2) using apply()
df_norm <- apply(df[,2:8], 2, function(x) log2(x/df$Signature.1))
# convert the result to a data frame and set column names
df_norm <- data.frame(df_norm)
colnames(df_norm) <- paste0("COSMIC.",colnames(df_norm))
df_norm$patientID=substr(row.names(df_norm),9,12)%>%as.integer()
df_norm=df_norm[,c("patientID","COSMIC.Signature.2","COSMIC.Signature.3","COSMIC.Signature.6","COSMIC.Signature.7","COSMIC.Signature.10","COSMIC.Signature.13")]
DNA=left_join(DNA,df_norm,by="patientID")
cols <- c("COSMIC.Signature.2","COSMIC.Signature.3","COSMIC.Signature.6",
          "COSMIC.Signature.7","COSMIC.Signature.10","COSMIC.Signature.13")
for (col in cols) {
  min_val <- min(DNA[[col]], na.rm = TRUE)   # 该列最小值
  DNA[[col]][is.na(DNA[[col]])] <- min_val  # 替换缺失值
}
colnames(DNA)[2:ncol(DNA)]=paste0("DNA_",colnames(DNA)[2:ncol(DNA)])

RNA=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
RNA_name=c("mRNA-ESR1","mRNA-ERBB2","mRNA-CD8A","mRNA-MKI67","mRNA-PGR","Apoptosis",                         
       "Lysosome","Endocytosis","Oxidative_phosphorylation","Purine_metabolism","Citrate_cycles",
       "Glutathione_metabolism","Fatty_acid_metabolism","Glycolysis","Hypoxia","EMT",                               
       "Exosome","sspbc.subtype","Taxane_response","HER2DX_IGG",                        
       "HER2DX_prolif","HER2DX_luminal","HER2DX_HER2_amplicon","HER2DX_pCR_likelihood_score","pik3ca_sig")
RNA <- RNA[, c("patientID",RNA_name)]
str(RNA)
RNA$patientID=as.double(RNA$patientID)
colnames(RNA)[2:ncol(RNA)]=paste0("RNA_",colnames(RNA)[2:ncol(RNA)])
ADC_traficking=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS9/ADC_traficking.rds")
ADC_traficking$RNA_ADC_trafficking=ADC_traficking$ADC_traficking
ADC_traficking$ADC_traficking=NULL
RNA=left_join(RNA,ADC_traficking,by="patientID")

#immune
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
data[,variable]=scale(data[,variable])
data=data[,c("patientID",variable)]
colnames(data)[2:ncol(data)]=paste0("RNA_",colnames(data)[2:ncol(data)])
RNA=left_join(RNA,data,by="patientID")
#DNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
variable=c("A01","meanHED","lohhla","TCRA.tcell.fraction.adj","Neoantigen_DNA")
data=data[,c("patientID",variable)]
data[,"lohhla"]=lapply(data[,"lohhla"], function(x) ifelse(x == 1, "TRUE", "FALSE"))
data[,c("A01")]=lapply(data[,c("A01")], function(x) ifelse(x == "Yes", "TRUE", "FALSE"))
colnames(data)[2]=paste0("HLA_Supertype_",colnames(data)[2])
colnames(data)[2:ncol(data)]=paste0("DNA_",colnames(data)[2:ncol(data)])
DNA=left_join(DNA,data,by="patientID")
DNA$patientID=as.double(DNA$patientID)
# image
image=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/image_metrics_PREDIX_HER2.rds")
colnames(image)[2:ncol(image)]=paste0("WSI_",colnames(image)[2:ncol(image)])
# protein
library(data.table);library(tidyverse);library(ggpubr)
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
MS=MS[c("RAB11FIP1","PPFIA1","CTTN","RPL19","MED1","CDK12","PPP1R1B",
        "SLC12A2","RAB11B","EEA1","ARL1","FLOT1","VAMP3","RAB5C"),]
MS=t(MS)%>%as.data.frame()
MS$patientID=row.names(MS)%>%as.double()
Prot=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
Prot=Prot[,c("patientID","ERBB2","GRB7","MIEN1","HER2_amplicon","ERBB2_PG")]
Prot=left_join(Prot,MS,by="patientID")
colnames(Prot)[2:ncol(Prot)]=paste0("Prot_",colnames(Prot)[2:ncol(Prot)])
#merge
data=left_join(clin,RNA,by='patientID')%>%left_join(DNA,by='patientID')%>%left_join(Prot,by='patientID')%>%left_join(image,by='patientID')%>%as.data.frame()
data$
write.table(data,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2_withNA.txt",quote = F,row.names =F,sep="\t")
#levels(data$RNA_sspbc.subtype) <- c(levels(data$RNA_sspbc.subtype), "Unknown")
#data$RNA_sspbc.subtype[is.na(data$RNA_sspbc.subtype)] <- "Unknown"
#data[is.na(data)]="Unknown"
no_na_data <- data[complete.cases(data), ]
nrow(no_na_data)    

write.table(no_na_data,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")
write.table(no_na_data[,colnames(clin)],file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")
write.table(no_na_data[,colnames(RNA)],file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/RNA_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")
write.table(no_na_data[,colnames(DNA)],file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/DNA_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")
write.table(no_na_data[,colnames(Prot)],file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/Prot_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")
write.table(no_na_data[,colnames(image)],file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/Image_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")

clin=data[,colnames(clin)]
clin <- clin[complete.cases(clin), ]
write.table(clin,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")

RNA=data[,colnames(RNA)]
RNA <- RNA[complete.cases(RNA), ]
write.table(RNA,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/RNA_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")

DNA=data[,colnames(DNA)]
DNA <- DNA[complete.cases(DNA), ]
write.table(DNA,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/DNA_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")

Prot=data[,colnames(Prot)]
Prot <- Prot[complete.cases(Prot), ]
write.table(Prot,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/Prot_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")

image=data[,colnames(image)]
image <- image[complete.cases(image), ]
write.table(image,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/Image_curated_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")


