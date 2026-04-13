############################################
#################BRCA driver################
############################################
library(data.table)
library(tidyverse)
library(maftools)
library(tableone)
library(readxl)
mutect=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_curated.txt")
maf=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_vaf0.05_curated.rds")
table(maf$Variant_Classification)
library(maftools)
laml = read.maf(maf = maf)
laml.mutload = tcgaCompare(maf = laml, cohortName = 'PREDIX HER2', logscale = TRUE, capture_size = 50)

vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                 "In_Frame_Ins", "Missense_Mutation")
freq=maf%>%group_by(Tumor_Sample_Barcode)%>%summarise(n = n())
median(freq$n)
maf=maf%>%filter(Variant_Classification%in%vc.nonSilent)
oncokb=maf%>%filter(ONCOGENIC%in%c("Oncogenic","Likely Oncogenic")) #"Likely Oncogenic"
oncokb$Hugo_Symbol%>%table()

HER_pathway=c("EGFR","ERBB2","ERBB3","ERBB4") # mut
PIK3_AKT_pathway=c("PIK3CA","PIK3R1","PTEN","AKT1") # PTEN del
MAPK_ERK_pathway=c("MAPK","MAP2K","MAP3K1","MAP2K4","NF1") # NF1 del
CDK_RB_pathway=c("CCND1","CDK4","CDK6","RB1") # CCND1 and CDK4/6 AMP  RB1 Del
# select most frequent mutation
gene_counts <- table(maf[["Hugo_Symbol"]])
gene_counts_df <- as.data.frame(gene_counts)
BRCA=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Resource/41586_2016_BFnature17676_MOESM47_ESM/nature17676-s3/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx",sheet=1)
BRCA=unique(BRCA$Gene)
# select gene for waterfall plot
BRCA=intersect(BRCA,gene_counts_df$Var1[gene_counts_df$Freq>9])
#BRCA=union(BRCA,gene_counts_df$Var1[gene_counts_df$Freq>10])
gene_counts_df=gene_counts_df[gene_counts_df$Var1%in%BRCA,]
saveRDS(BRCA,file="E:/Projects/PREDIX_HER2/Multimodal/Resource/predix_her2_freqORdriver.rds")
#BRCA = fread("E:/Projects/Collaboration/BEVPAC/CUTseq/genelist_nik-zainal-etal.tsv") 
df=maf%>%filter(Hugo_Symbol%in%c(BRCA,HER_pathway,PIK3_AKT_pathway,MAPK_ERK_pathway,CDK_RB_pathway))
table(df$Hugo_Symbol)
mutation=df[,c("Hugo_Symbol","Tumor_Sample_Barcode")]
mutation$value=1
mutation$Tumor_Sample_Barcode=as.factor(mutation$Tumor_Sample_Barcode)
data = reshape2::dcast(mutation, Tumor_Sample_Barcode~Hugo_Symbol,value.var = "value")
data[data==0]="WT"
data[data==1]="Mut"
data[data==2]="Mut"
data[data==3]="Mut"
data[data==4]="Mut"
data[data==5]="Mut"
data[data==6]="Mut"
colnames(data)[1]="sampleID"
colnames(data)[2:ncol(data)]=paste0("coding_mutation_",colnames(data)[2:ncol(data)])
driver=data.frame(sampleID=unique(mutect$Tumor_Sample_Barcode))
driver=left_join(driver,data,by="sampleID")
driver[is.na(driver)]=0
driver[driver=='WT']=0
driver[driver=='Mut']=1
driver$patientID=substr(driver$sampleID,9,12)
####For pathway only oncokb annotated###
# CUTseq gistic2 and WES as complementary
CNA=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/CUTseq_gene_baseline_complemental.rds")
CNA$patientID=row.names(CNA)
# HER family pathway #
driver$coding_mutation_HER_pathway=0
driver$coding_mutation_HER_pathway[driver$sampleID%in%maf$Tumor_Sample_Barcode[maf$Hugo_Symbol%in%c("EGFR","ERBB2","ERBB3","ERBB4")]]=1
# PIK3_AKT_pathway
PIK3CA=CNA$patientID[CNA$PIK3CA%in%c(2)] # AMPLIFICATION
AKT1=CNA$patientID[CNA$AKT1%in%c(2)] # AMPLIFICATION
PTEN=CNA$patientID[CNA$PTEN%in%c(-2)] # deletion
driver$coding_mutation_PIK3_AKT_pathway=0
driver$coding_mutation_PIK3_AKT_pathway[driver$sampleID%in%maf$Tumor_Sample_Barcode[maf$Hugo_Symbol%in%c("PIK3CA","PIK3R1","PTEN","AKT1")]]=1
driver$coding_mutation_PIK3_AKT_pathway[driver$patientID%in%c(PIK3CA,AKT1,PTEN)]=1
# MAPK_ERK_pathway
NF1=CNA$patientID[CNA$NF1%in%c(-2)] # deletion
driver$coding_mutation_MAPK_ERK_pathway=0
driver$coding_mutation_MAPK_ERK_pathway[driver$sampleID%in%maf$Tumor_Sample_Barcode[maf$Hugo_Symbol%in%c("MAPK","MAP2K","MAP3K1","MAP2K4","NF1")]]=1
driver$coding_mutation_MAPK_ERK_pathway[driver$patientID%in%c(NF1)]=1
table(driver$coding_mutation_MAPK_ERK_pathway) 
# CDK_RB_pathway
CCND1=CNA$patientID[CNA$CCND1%in%c(2)] # AMPLIFICATION
CDK4=CNA$patientID[CNA$CDK4%in%c(2)] # AMPLIFICATION
CDK6=CNA$patientID[CNA$CDK6%in%c(2)] # AMPLIFICATION
RB1=CNA$patientID[CNA$RB1%in%c(-2)] # deletion
driver$coding_mutation_CDK_RB_pathway=0
driver$coding_mutation_CDK_RB_pathway[driver$sampleID%in%maf$Tumor_Sample_Barcode[maf$Hugo_Symbol%in%c("CCND1","CDK4","CDK6","RB1")]]=1
driver$coding_mutation_CDK_RB_pathway[driver$patientID%in%c(CCND1,CDK4,CDK6,RB1)]=1
table(driver$coding_mutation_CDK_RB_pathway)
#TP53_oncokb
driver$coding_mutation_TP53_oncokb=0
driver$coding_mutation_TP53_oncokb[driver$sampleID%in%oncokb$Tumor_Sample_Barcode[oncokb$Hugo_Symbol%in%c("TP53")]]=1
#PIK3CA_oncokb
driver$coding_mutation_PIK3CA_oncokb=0
driver$coding_mutation_PIK3CA_oncokb[driver$sampleID%in%oncokb$Tumor_Sample_Barcode[oncokb$Hugo_Symbol%in%c("PIK3CA")]]=1
#ERBB2_oncokb
driver$coding_mutation_ERBB2_oncokb=0
driver$coding_mutation_ERBB2_oncokb[driver$sampleID%in%oncokb$Tumor_Sample_Barcode[oncokb$Hugo_Symbol%in%c("ERBB2")]]=1
#GATA3_oncokb
driver$coding_mutation_GATA3_oncokb=0
driver$coding_mutation_GATA3_oncokb[driver$sampleID%in%oncokb$Tumor_Sample_Barcode[oncokb$Hugo_Symbol%in%c("GATA3")]]=1
# merge
BRCA=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Resource/predix_her2_freqORdriver.rds")
BRCA=paste0("coding_mutation_",BRCA)
driver=driver[,c("sampleID",BRCA,"coding_mutation_HER_pathway","coding_mutation_PIK3_AKT_pathway",
               "coding_mutation_MAPK_ERK_pathway","coding_mutation_CDK_RB_pathway","coding_mutation_TP53_oncokb","coding_mutation_PIK3CA_oncokb",
               "coding_mutation_GATA3_oncokb","coding_mutation_ERBB2_oncokb")]
driver[,2:ncol(driver)]=lapply(as.data.frame(driver[,2:ncol(driver)]),function(x) as.factor(x))
str(driver)
table(driver$coding_mutation_TP53_oncokb)
############################################
#################BRCA CNA CUTseq############
############################################
CNA_brca=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/CUTseq_gene_baseline_GISTIC_score_complemental.rds")
drivers=fread("E:/Projects/Collaboration/BEVPAC/CUTseq/genelist_nik-zainal-etal.tsv") 
gene=c("ERBB2","BRCA2","NCOR1","PIK3CA","RAB11FIP1","FADD","PPFIA1","CTTN","RPL19","MED1","CDK12","PPP1R1B","MIEN1","GRB7")

CNA_brca=CNA_brca[,gene]
colnames(CNA_brca)=paste0(colnames(CNA_brca),"_CNA")
CNA_brca$patientID=row.names(CNA_brca)
############################################
##########Tumor mutational signature########
############################################
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
offset=0.001
mutect=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_curated.txt")
sample=data.frame(sampleID=unique(mutect$Tumor_Sample_Barcode))
df$sampleID=row.names(df)
df=left_join(sample,df,by='sampleID')
df[df==0|is.na(df)]=offset
row.names(df)=df$sampleID
df$sampleID=NULL
# normalize columns using the first column (log2) using apply()
df_norm <- apply(df[,2:8], 2, function(x) log2(x/df$Signature.1))
# convert the result to a data frame and set column names
df_norm <- data.frame(df_norm)
colnames(df_norm) <- paste0("COSMIC.",colnames(df_norm))
df_norm$sampleID=row.names(df_norm)
#############################################
################PureCN CNA###################
#############################################
library(tidyverse);library(data.table);library(readxl)
loh=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_LOH_curated.csv")) 
table(loh$type,loh$M.flagged)
CNV=loh%>%filter(C!=2,M!=1)%>%filter(chr!="chrX")
# genome length https://www.ncbi.nlm.nih.gov/grc/human/data
length=248956422+242193529+198295559+190214555+181538259+170805979+159345973+145138636+138394717+133797422+135086622+
  133275309+114364328+107043718+101991189+90338345+83257441+80373285+58617616+64444167+46709983+50818468
CNV$genome_fraction=(CNV$end-CNV$start)/length
CNV=CNV%>%group_by(sampleID)%>%summarise(CNV_burden=sum(genome_fraction))
LOH=loh%>%filter(M.flagged!='TRUE')%>%filter(chr!="chrX") #flag indicating that M is unreliable
LOH=LOH%>%filter(type%in%c("LOH","WHOLE ARM LOH","COPY-NEUTRAL LOH","WHOLE ARM COPY-NEUTRAL LOH"))
LOH$genome_fraction=(LOH$end-LOH$start)/(1000000*3100)
LOH$type%>%table()
LOH_Del=LOH%>%filter(C%in%c("0","1"))%>%group_by(sampleID)%>%summarise(LOH_Del_burden=sum (genome_fraction))
LOH=data.frame(sampleID=unique(loh$sampleID))
CNA=left_join(LOH,LOH_Del,by="sampleID")%>%left_join(CNV,by="sampleID")
CNA[is.na(CNA)]=0
###########################################
############Subclone percentage############
###########################################
# From PureCN purity, copy number and SNV
mutect=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_vaf0.05_curated.rds")
maf=mutect
maf$vaf=maf$t_alt_count/maf$t_depth
somatic.snv=as.data.frame(maf)
# PureCN
purity=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity_curated.csv"))
colnames(purity)[2:3]=c("purity",'ploidy')
purity=purity[,c('sampleID','ploidy',"purity")]
purity=purity[purity$sampleID%in%unique(mutect$Tumor_Sample_Barcode),]
allele.cnv=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_LOH_curated.csv"))
allele.cnv=allele.cnv[!is.na(allele.cnv$type),]
allele.cnv$nMajor=allele.cnv$C-allele.cnv$M
allele.cnv$nMinor=allele.cnv$M
allele.cnv$chrom=allele.cnv$chr
allele.cnv=allele.cnv[,c("sampleID",'chrom',"start","end","nMajor","nMinor")]
head(allele.cnv)
allele.cnv=allele.cnv[allele.cnv$sampleID%in%unique(mutect$Tumor_Sample_Barcode),]
save(somatic.snv,purity,allele.cnv,file ="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PREDIX_HER2_Baseline_Clonal_analyses.RData") 
library(foreach)
#############################################
# Function to calculate cancer cell fraction
# Adapted from Mcgranahan et al, Science
absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number)
{
  f.function <- function (c,purity,local.copy.number)
  {
    return(min(c((purity*c) / (2*(1-purity) + purity*local.copy.number),1)))
  }
  x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number))
  if(is.na(min(x))){
    d = data.frame(
      ccf.lower.ci = NA, ccf.est = NA, ccf.upper.ci = NA,
      prob.subclonal=NA,prob.clonal=NA, 
      is.clonal=NA,is.clonal.byprob=NA
    )
    return (d)
  }
  if(min(x)==0)
  {
    x[length(x)] <- 1
  }
  names(x)       <- seq(0.01,1,length.out=100)
  sub.cint <- function(x, prob = 0.95,n.alt,depth) {
    xnorm   <- x/sum(x)
    xsort   <- sort(xnorm, decreasing = TRUE)
    xcumLik <- cumsum(xsort)
    n = sum(xcumLik < prob) + 1
    LikThresh <- xsort[n]
    cint  <- x[xnorm >= LikThresh]
    all   <- as.numeric(names(x))
    cellu <- as.numeric(names(cint))
    l.t   <- cellu[1]
    r.t   <- cellu[length(cellu)]
    m     <- cellu[which.max(cint)]
    
    prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative=‘less’)$p.val
    prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative=‘greater’)$p.val
    
    data.frame(
      ccf.lower.ci = l.t, ccf.est = m, ccf.upper.ci = r.t,
      prob.subclonal=prob.subclonal,prob.clonal=prob.clonal, 
      # Define clonality strictly as a confidence interval overlapping 1
      is.clonal=(l.t <= 1 & r.t >= 1),
      # Define clonality differently as prob.clonal > prob.subclonal (only for comparison purposes)
      is.clonal.byprob=(prob.clonal > prob.subclonal)
    )
  }
  return(sub.cint(x,n.alt=n.alt,depth=depth))
}
#
load("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PREDIX_HER2_Baseline_Clonal_analyses.RData")
somatic.snv=somatic.snv%>%filter(Variant_Type=="SNP")
unique(purity$sampleID);unique(allele.cnv$sampleID);unique(somatic.snv$Tumor_Sample_Barcode)
n = nrow(purity)
clonal_results = foreach (i=1:n, pat=purity$sampleID,.combine=rbind) %do% {
  cat(pat, "...\n")
  ## snv information
  snv=subset(somatic.snv,Tumor_Sample_Barcode==pat)
  rownames(snv) = paste(snv$Chromosome, snv$Start_Position, sep="_")
  snv. = data.frame(
    chrom = snv$Chromosome,
    position = snv$Start_Position,
    readcount = snv$t_ref_count + snv$t_alt_count,
    alt_count = snv$t_alt_count,
    vaf = snv$t_alt_count/(snv$t_ref_count + snv$t_alt_count)
  )
  if(nrow(snv.)==0){
    return (NULL)
  }
  if(!grepl("chr",snv$Chromosome[1])){
    snv.$chrom = paste("chr", snv$Chromosome,sep="")
  }
  ## cnv information
  cnv = subset(allele.cnv, sampleID==pat)
  cnv$nTotal = cnv$nMajor+cnv$nMinor
  
  indices = sapply(1:nrow(snv.), function(x) 
    which(cnv$chrom==snv.$chrom[x] & cnv$start<=snv.$position[x] & cnv$end>=snv.$position[x])[1])
  snv.$cnTotal = cnv[indices, "nTotal"]
  snv.$cnMajor = cnv[indices, "nMajor"]
  snv.$cnMinor = cnv[indices, "nMinor"]
  
  rho = purity$purity[i] ## tumor purity
  snv.$nmut = snv.$vaf * (rho * snv.$cnTotal + (1-rho) * 2) / rho ## multiplicity
  ## calculate CCF
  ccf = lapply(1:nrow(snv.), function (x){
    absolute.cancer.cell.fraction(
      n.alt = snv.$alt_count[x],
      depth = snv.$readcount[x],
      purity = rho,
      local.copy.number=snv.$cnTotal[x]
    )
  })
  ccf <- data.table::rbindlist(ccf)
  return(cbind(snv[,1:(ncol(snv)-4)], snv.[,3:ncol(snv.)],ccf))
  
}
saveRDS(clonal_results,file="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PREDIX_HER2_Baseline_Clonal_analyses.rds")
## final results
head(clonal_results)
table(clonal_results$is.clonal)
freq=clonal_results%>%group_by(Tumor_Sample_Barcode)%>%summarise(n = n())
subclone=data.frame(sampleID=unique(clonal_results[,"Tumor_Sample_Barcode"]))  
subclone$subclone_per=NA
for (x in subclone$sampleID) {
  sub=clonal_results[clonal_results$Tumor_Sample_Barcode==x,]
  subclone_per=nrow(sub[sub$is.clonal=="FALSE",])/nrow(sub)
  subclone[subclone$sampleID==x,"subclone_per"]=subclone_per
}

subclone
ggplot(subclone, aes(x=subclone_per)) + geom_histogram()
subclone
############################################
##########Tumor mutational burden###########
############################################
# Uniform TMB estimation:  coding, non-synonymous single-nucleotide polymorphisms and insertions and deletions
# VAF >5%  genome size 36.8 Mb
library(data.table)
library(tidyverse)
library(maftools)
library(tableone)
clonal_results=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PREDIX_HER2_Baseline_Clonal_analyses.rds")
loh=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_gene_curated.csv")) 
loh$Tumor_Sample_Barcode=loh$sampleID;loh$Hugo_Symbol=loh$gene.symbol
clonal_results=left_join(clonal_results,loh,by=c("Tumor_Sample_Barcode","Hugo_Symbol"))
vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                 "In_Frame_Ins", "Missense_Mutation")
freq=clonal_results%>%dplyr::group_by(Tumor_Sample_Barcode)%>%dplyr::summarise(n = n())%>%as.data.frame()
total=data.frame(sampleID=freq$Tumor_Sample_Barcode,totalTMB=freq$n/41.2)
maf=clonal_results%>%filter(Variant_Classification%in%vc.nonSilent,vaf>=0.05,n_depth>25,t_depth>25,is.na(gnomAD_AF)|gnomAD_AF<0.01)
freq=maf%>%dplyr::group_by(Tumor_Sample_Barcode)%>%dplyr::summarise(n = n())%>%as.data.frame()
tmb=data.frame(sampleID=freq$Tumor_Sample_Barcode,TMB_uniform=freq$n/41.2)
maf_clone=maf%>%filter(is.clonal=="TRUE")
maf_subclone=maf%>%filter(is.clonal=="FALSE")
freq=maf_clone%>%group_by(Tumor_Sample_Barcode)%>%summarise(n = n())%>%as.data.frame()
freq_sub=maf_subclone%>%group_by(Tumor_Sample_Barcode)%>%summarise(n = n())%>%as.data.frame()
tmb_clone=data.frame(sampleID=freq$Tumor_Sample_Barcode,TMB_clone=freq$n/41.2)
tmb_subclone=data.frame(sampleID=freq_sub$Tumor_Sample_Barcode,TMB_subclone=freq_sub$n/41.2)

tmb=left_join(total,tmb,by="sampleID")%>%left_join(tmb_clone,by="sampleID")%>%
  left_join(tmb_subclone,by="sampleID")
tmb[is.na(tmb)]=0
###########################################
###############HRD score###################
###########################################
library("scarHRD");library(foreach) # read in allele specific CN
# a=scar_score("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/test2.txt",reference = "grch38",seqz=FALSE)
# fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/test2.txt")
CN=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_LOH_curated.csv"),stringsAsFactors=FALSE)
unique(CN$sampleID)
CN=CN%>%filter(num.snps!=0) 
CN$C=as.numeric(CN$C);CN$M=as.numeric(CN$M)
CN=data.frame(SampleID=CN$sampleID,
              Chromosome=CN$chr, Start_position=CN$start,End_position=CN$end,
              total_cn=CN$C,
              A_cn=CN$M,B_cn=CN$C-CN$M,ploidy=NA) # A_cn minor B_cn major
sampleID=unique(CN$SampleID)
HRD= foreach (i=sampleID,.combine=rbind) %do% {
  cat(i, "...\n")  
  cn=filter(CN,SampleID==i)
  write.table(cn,file="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/scarHRD.input.txt",quote = F,row.names =T,sep = "\t")
  HRD=as.data.frame(scar_score("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/scarHRD.input.txt",reference = "grch38",seqz=FALSE,ploidy =NULL)) 
}
HRD$sampleID=sampleID
str(HRD)
HRD[,c("HRD","Telomeric AI","LST")]=NULL
colnames(HRD)=c("HRD","sampleID")
######################################################
##################Data Integration####################
######################################################
driver$patientID=substr(driver$sampleID,9,12)
genomic=left_join(driver,CNA_brca,by="patientID")%>%
  left_join(df_norm,by="sampleID")%>%
  left_join(CNA,by="sampleID")%>%
  left_join(subclone,by="sampleID")%>%
  left_join(tmb,by="sampleID")%>%
  left_join(HRD,by="sampleID")

cols <- c("subclone_per","totalTMB","TMB_uniform","TMB_clone")
for (col in cols) {
  genomic[[col]][is.na(genomic[[col]])] <- 0
}

cols <- c("COSMIC.Signature.2","COSMIC.Signature.3","COSMIC.Signature.6","COSMIC.Signature.7","COSMIC.Signature.10","COSMIC.Signature.13","COSMIC.Signature.15")
genomic[genomic$totalTMB*41.2<10,cols]=NA 
#genomic[genomic$totalTMB*41.2<10,"TMB_clone"]=NA 

str(genomic)
colnames(genomic)
#genomic$patientID=substr(genomic$sampleID,9,12)%>%as.double()
genomic<- genomic[, c("sampleID", "patientID", setdiff(names(genomic), c("sampleID", "patientID")))]
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=genomic[genomic$patientID%in%clin$patientID,]
genomic$patientID=as.integer(genomic$patientID)
saveRDS(genomic,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
write.table(genomic,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt",quote = F,row.names =F,sep="\t")












