# Title     : TODO
# Objective : TODO
# Created by: Chul
# Created on: 2020-09-09

# GLOBAL_var
wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"

nLinkType = "Full"

# Load libraries
#install.packages("viridis")  # Install
library("viridis")
library("OmicCircos")
options(stringsAsFactors=FALSE)
library("stringr")



#oNAME = "chrX.pdf"
setwd(wPATH)
iPATH             = paste0(wPATH,"BUSCO_circos_input_byChr/")
oPATH = paste0(wPATH,"output/BUSCO_circos_output_byChr/")
ChrID_list = c()
flist = list.files(iPATH)
for(fNAME in flist){
if(str_sub(fNAME, 1,9)==paste0(nLinkType,"Link_")){
fNAME = strsplit(fNAME,"_")[[1]][2]
ChrID = strsplit(fNAME,".tsv")[[1]][1]
ChrID_list = c(ChrID_list,ChrID)
}
}
color_target_sID_chr = "#00c4a7"
color_others_sID_chr = "#E7B800"
iCNT_sID = 25 #pal_jco()(iCNT_sID)
#color_code_chr = c("black",rainbow(iCNT_sID))
#color_code_link = rainbow(iCNT_sID,alpha = .2)
name = "Zissou1"
color_code_chr = c("black",plasma(iCNT_sID))
color_code_link = plasma(iCNT_sID, alpha=.1)


# Function
circos_plot <- function(assemblyID,fNAME_chr_size,fNAME_chr_name,fNAME_link,fNAME_color_chr,fNAME_color_link,oNAME_pdf) 
{
# Load data
chr_size_df <- read.table(fNAME_chr_size,sep="\t",header=T)
chr_names_df <- read.table(fNAME_chr_name,sep="\t",header=T)
link_df <- read.table(fNAME_link,sep="\t",header=T)
color_chr_df <- read.table(fNAME_color_chr,sep="\t",header=T)
color_link_df <- read.table(fNAME_color_link,sep="\t",header=T)
# Data pre-processing
chr_order <- chr_size_df$chrom
seg.name <- chr_order
chr_db <- segAnglePo(chr_size_df, seg=seg.name)
color_chr <- color_code_chr[color_chr_df$color]
color_link <- color_code_link[color_link_df$color] 
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main=assemblyID);
circos(cir=chr_db,R=180,W=10,type="chr", col=color_chr, print.chr.lab=FALSE,scale=FALSE,order=chr_order,lwd=.5,cex=0.1)
circos(cir=chr_db,R=190,W=0,mapping=chr_names_df,type="label", side="out",col=color_chr,order=chr_order, cex=0.47, lwd = 0)
circos(cir=chr_db,R=170,W=0,type="link", mapping=link_df, col = color_link, lwd=.01)
}
BUSCO_circos_plot <- function(ChrID){
fNAME_chr_size = paste0(iPATH,nLinkType,"ChrSize_",ChrID,".tsv",collapse = NULL)
fNAME_chr_name = paste0(iPATH,nLinkType,"ChrName_",ChrID,".tsv",collapse = NULL)
fNAME_color_chr = paste0(iPATH,nLinkType,"ChrColor_",ChrID,".tsv",collapse = NULL)
fNAME_link = paste0(iPATH,nLinkType,"Link_",ChrID,".tsv",collapse = NULL)
fNAME_color_link = paste0(iPATH,nLinkType,"LinkColor_",ChrID,".tsv",collapse = NULL)
oNAME_pdf = paste0(oPATH,"BUSCO_circos_plot_byChr_",nLinkType,"_",ChrID,".pdf",collapse = NULL)
pdf(oNAME_pdf,6,6)
#pdf(oNAME_pdf,10,10)
par(mar=c(1.5,1.5,1.5,1.5))
circos_plot(ChrID,fNAME_chr_size,fNAME_chr_name,fNAME_link,fNAME_color_chr,fNAME_color_link,oNAME_pdf) 
dev.off()
}

#pdf(oNAME,10,10)
nLinkType_vec = c("Full","Core")
for (nLinkType in nLinkType_vec){
for (ChrID in ChrID_list)
{
BUSCO_circos_plot(ChrID)
}
}

#BUSCO_circos_plot("chr1")
