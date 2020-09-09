RefID = "hg38"
iCNT_RefChr = 35
wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"
col_vec_comp = rainbow(iCNT_RefChr, alpha = 0.2)
col_vec_dnaseg = rainbow(iCNT_RefChr, alpha = 1) 

library(genoPlotR)
setwd(wPATH)

iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID



                                           
# function
# generate_dnaseg
make_dnaseq_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/dnaseg/dnaseg_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
if(sID==RefID){
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=col_vec_dnaseg[data$fill])
}else{
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=data$fill)
}
dna_segX <- dna_seg(df,gene_type = "blocks")
return(dna_segX)
}

# generate_comp
make_comp_eachID <- function(compNumb){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/comparison/comparison",compNumb,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(start1=data$start1, end1=data$end1, start2=data$start2, end2 = data$end2, col = col_vec_comp[data$col],direction=1)
comparisonX = comparison(df)
return(comparisonX)
}

dna_segs = list(
make_dnaseq_eachID("hg38"),
make_dnaseq_eachID("aRhiBiv1"),
make_dnaseq_eachID("bCalAnn1"),
make_dnaseq_eachID("bStrHab1"),
make_dnaseq_eachID("bTaeGut2"),
make_dnaseq_eachID("fAnaTes1"),
make_dnaseq_eachID("fArcCen1"),
make_dnaseq_eachID("fAstCal1"),
make_dnaseq_eachID("fCotGob3"),
make_dnaseq_eachID("fGouWil2"),
make_dnaseq_eachID("fMasArm1"),
make_dnaseq_eachID("mLynCan4"),
make_dnaseq_eachID("mOrnAna1"),
make_dnaseq_eachID("mPhyDis1"),
make_dnaseq_eachID("mRhiFer1"),
make_dnaseq_eachID("rGopEvg1"),
make_dnaseq_eachID("sAmbRad1")
)
names(dna_segs)=sID_list

comparisons = list(
make_comp_eachID("01"),
make_comp_eachID("02"),
make_comp_eachID("03"),
make_comp_eachID("04"),
make_comp_eachID("05"),
make_comp_eachID("06"),
make_comp_eachID("07"),
make_comp_eachID("08"),
make_comp_eachID("09"),
make_comp_eachID("10"),
make_comp_eachID("11"),
make_comp_eachID("12"),
make_comp_eachID("13"),
make_comp_eachID("14"),
make_comp_eachID("15"),
make_comp_eachID("16")
)

pdf(paste0(wPATH,"output/BUSCO_genoPlotR_output/BUSCO_genoPlotR_colorred_by_",RefID,".pdf"),6,9)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,main=paste0("CoLink_Colorred_by_",RefID))
dev.off()

RefID = "mOrnAna1"
iCNT_RefChr = 43
wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"
col_vec_comp = rainbow(iCNT_RefChr, alpha = 0.2)
col_vec_dnaseg = rainbow(iCNT_RefChr, alpha = 1) 

library(genoPlotR)
setwd(wPATH)

iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID



                                           
# function
# generate_dnaseg
make_dnaseq_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/dnaseg/dnaseg_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
if(sID==RefID){
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=col_vec_dnaseg[data$fill])
}else{
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=data$fill)
}
dna_segX <- dna_seg(df,gene_type = "blocks")
return(dna_segX)
}

# generate_comp
make_comp_eachID <- function(compNumb){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/comparison/comparison",compNumb,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(start1=data$start1, end1=data$end1, start2=data$start2, end2 = data$end2, col = col_vec_comp[data$col],direction=1)
comparisonX = comparison(df)
return(comparisonX)
}

dna_segs = list(
make_dnaseq_eachID("hg38"),
make_dnaseq_eachID("aRhiBiv1"),
make_dnaseq_eachID("bCalAnn1"),
make_dnaseq_eachID("bStrHab1"),
make_dnaseq_eachID("bTaeGut2"),
make_dnaseq_eachID("fAnaTes1"),
make_dnaseq_eachID("fArcCen1"),
make_dnaseq_eachID("fAstCal1"),
make_dnaseq_eachID("fCotGob3"),
make_dnaseq_eachID("fGouWil2"),
make_dnaseq_eachID("fMasArm1"),
make_dnaseq_eachID("mLynCan4"),
make_dnaseq_eachID("mOrnAna1"),
make_dnaseq_eachID("mPhyDis1"),
make_dnaseq_eachID("mRhiFer1"),
make_dnaseq_eachID("rGopEvg1"),
make_dnaseq_eachID("sAmbRad1")
)
names(dna_segs)=sID_list

comparisons = list(
make_comp_eachID("01"),
make_comp_eachID("02"),
make_comp_eachID("03"),
make_comp_eachID("04"),
make_comp_eachID("05"),
make_comp_eachID("06"),
make_comp_eachID("07"),
make_comp_eachID("08"),
make_comp_eachID("09"),
make_comp_eachID("10"),
make_comp_eachID("11"),
make_comp_eachID("12"),
make_comp_eachID("13"),
make_comp_eachID("14"),
make_comp_eachID("15"),
make_comp_eachID("16")
)

pdf(paste0(wPATH,"output/BUSCO_genoPlotR_output/BUSCO_genoPlotR_colorred_by_",RefID,".pdf"),6,9)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,main=paste0("CoLink_Colorred_by_",RefID))
dev.off()

RefID = "bTaeGut2"
iCNT_RefChr = 44
wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"
col_vec_comp = rainbow(iCNT_RefChr, alpha = 0.2)
col_vec_dnaseg = rainbow(iCNT_RefChr, alpha = 1) 

library(genoPlotR)
setwd(wPATH)

iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID



                                           
# function
# generate_dnaseg
make_dnaseq_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/dnaseg/dnaseg_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
if(sID==RefID){
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=col_vec_dnaseg[data$fill])
}else{
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=data$fill)
}
dna_segX <- dna_seg(df,gene_type = "blocks")
return(dna_segX)
}

# generate_comp
make_comp_eachID <- function(compNumb){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/comparison/comparison",compNumb,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(start1=data$start1, end1=data$end1, start2=data$start2, end2 = data$end2, col = col_vec_comp[data$col],direction=1)
comparisonX = comparison(df)
return(comparisonX)
}

dna_segs = list(
make_dnaseq_eachID("hg38"),
make_dnaseq_eachID("aRhiBiv1"),
make_dnaseq_eachID("bCalAnn1"),
make_dnaseq_eachID("bStrHab1"),
make_dnaseq_eachID("bTaeGut2"),
make_dnaseq_eachID("fAnaTes1"),
make_dnaseq_eachID("fArcCen1"),
make_dnaseq_eachID("fAstCal1"),
make_dnaseq_eachID("fCotGob3"),
make_dnaseq_eachID("fGouWil2"),
make_dnaseq_eachID("fMasArm1"),
make_dnaseq_eachID("mLynCan4"),
make_dnaseq_eachID("mOrnAna1"),
make_dnaseq_eachID("mPhyDis1"),
make_dnaseq_eachID("mRhiFer1"),
make_dnaseq_eachID("rGopEvg1"),
make_dnaseq_eachID("sAmbRad1")
)
names(dna_segs)=sID_list

comparisons = list(
make_comp_eachID("01"),
make_comp_eachID("02"),
make_comp_eachID("03"),
make_comp_eachID("04"),
make_comp_eachID("05"),
make_comp_eachID("06"),
make_comp_eachID("07"),
make_comp_eachID("08"),
make_comp_eachID("09"),
make_comp_eachID("10"),
make_comp_eachID("11"),
make_comp_eachID("12"),
make_comp_eachID("13"),
make_comp_eachID("14"),
make_comp_eachID("15"),
make_comp_eachID("16")
)

pdf(paste0(wPATH,"output/BUSCO_genoPlotR_output/BUSCO_genoPlotR_colorred_by_",RefID,".pdf"),6,9)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,main=paste0("CoLink_Colorred_by_",RefID))
dev.off()

RefID = "rGopEvg1"
iCNT_RefChr = 35
wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"
col_vec_comp = rainbow(iCNT_RefChr, alpha = 0.2)
col_vec_dnaseg = rainbow(iCNT_RefChr, alpha = 1) 

library(genoPlotR)
setwd(wPATH)

iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID



                                           
# function
# generate_dnaseg
make_dnaseq_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/dnaseg/dnaseg_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
if(sID==RefID){
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=col_vec_dnaseg[data$fill])
}else{
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=data$fill)
}
dna_segX <- dna_seg(df,gene_type = "blocks")
return(dna_segX)
}

# generate_comp
make_comp_eachID <- function(compNumb){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/comparison/comparison",compNumb,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(start1=data$start1, end1=data$end1, start2=data$start2, end2 = data$end2, col = col_vec_comp[data$col],direction=1)
comparisonX = comparison(df)
return(comparisonX)
}

dna_segs = list(
make_dnaseq_eachID("hg38"),
make_dnaseq_eachID("aRhiBiv1"),
make_dnaseq_eachID("bCalAnn1"),
make_dnaseq_eachID("bStrHab1"),
make_dnaseq_eachID("bTaeGut2"),
make_dnaseq_eachID("fAnaTes1"),
make_dnaseq_eachID("fArcCen1"),
make_dnaseq_eachID("fAstCal1"),
make_dnaseq_eachID("fCotGob3"),
make_dnaseq_eachID("fGouWil2"),
make_dnaseq_eachID("fMasArm1"),
make_dnaseq_eachID("mLynCan4"),
make_dnaseq_eachID("mOrnAna1"),
make_dnaseq_eachID("mPhyDis1"),
make_dnaseq_eachID("mRhiFer1"),
make_dnaseq_eachID("rGopEvg1"),
make_dnaseq_eachID("sAmbRad1")
)
names(dna_segs)=sID_list

comparisons = list(
make_comp_eachID("01"),
make_comp_eachID("02"),
make_comp_eachID("03"),
make_comp_eachID("04"),
make_comp_eachID("05"),
make_comp_eachID("06"),
make_comp_eachID("07"),
make_comp_eachID("08"),
make_comp_eachID("09"),
make_comp_eachID("10"),
make_comp_eachID("11"),
make_comp_eachID("12"),
make_comp_eachID("13"),
make_comp_eachID("14"),
make_comp_eachID("15"),
make_comp_eachID("16")
)

pdf(paste0(wPATH,"output/BUSCO_genoPlotR_output/BUSCO_genoPlotR_colorred_by_",RefID,".pdf"),6,9)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,main=paste0("CoLink_Colorred_by_",RefID))
dev.off()

RefID = "aRhiBiv1"
iCNT_RefChr = 28
wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"
col_vec_comp = rainbow(iCNT_RefChr, alpha = 0.2)
col_vec_dnaseg = rainbow(iCNT_RefChr, alpha = 1) 

library(genoPlotR)
setwd(wPATH)

iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID



                                           
# function
# generate_dnaseg
make_dnaseq_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/dnaseg/dnaseg_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
if(sID==RefID){
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=col_vec_dnaseg[data$fill])
}else{
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=data$fill)
}
dna_segX <- dna_seg(df,gene_type = "blocks")
return(dna_segX)
}

# generate_comp
make_comp_eachID <- function(compNumb){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/comparison/comparison",compNumb,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(start1=data$start1, end1=data$end1, start2=data$start2, end2 = data$end2, col = col_vec_comp[data$col],direction=1)
comparisonX = comparison(df)
return(comparisonX)
}

dna_segs = list(
make_dnaseq_eachID("hg38"),
make_dnaseq_eachID("aRhiBiv1"),
make_dnaseq_eachID("bCalAnn1"),
make_dnaseq_eachID("bStrHab1"),
make_dnaseq_eachID("bTaeGut2"),
make_dnaseq_eachID("fAnaTes1"),
make_dnaseq_eachID("fArcCen1"),
make_dnaseq_eachID("fAstCal1"),
make_dnaseq_eachID("fCotGob3"),
make_dnaseq_eachID("fGouWil2"),
make_dnaseq_eachID("fMasArm1"),
make_dnaseq_eachID("mLynCan4"),
make_dnaseq_eachID("mOrnAna1"),
make_dnaseq_eachID("mPhyDis1"),
make_dnaseq_eachID("mRhiFer1"),
make_dnaseq_eachID("rGopEvg1"),
make_dnaseq_eachID("sAmbRad1")
)
names(dna_segs)=sID_list

comparisons = list(
make_comp_eachID("01"),
make_comp_eachID("02"),
make_comp_eachID("03"),
make_comp_eachID("04"),
make_comp_eachID("05"),
make_comp_eachID("06"),
make_comp_eachID("07"),
make_comp_eachID("08"),
make_comp_eachID("09"),
make_comp_eachID("10"),
make_comp_eachID("11"),
make_comp_eachID("12"),
make_comp_eachID("13"),
make_comp_eachID("14"),
make_comp_eachID("15"),
make_comp_eachID("16")
)

pdf(paste0(wPATH,"output/BUSCO_genoPlotR_output/BUSCO_genoPlotR_colorred_by_",RefID,".pdf"),6,9)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,main=paste0("CoLink_Colorred_by_",RefID))
dev.off()

RefID = "sAmbRad1"
iCNT_RefChr = 70
wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"
col_vec_comp = rainbow(iCNT_RefChr, alpha = 0.2)
col_vec_dnaseg = rainbow(iCNT_RefChr, alpha = 1) 

library(genoPlotR)
setwd(wPATH)

iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID



                                           
# function
# generate_dnaseg
make_dnaseq_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/dnaseg/dnaseg_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
if(sID==RefID){
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=col_vec_dnaseg[data$fill])
}else{
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=data$fill)
}
dna_segX <- dna_seg(df,gene_type = "blocks")
return(dna_segX)
}

# generate_comp
make_comp_eachID <- function(compNumb){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/comparison/comparison",compNumb,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(start1=data$start1, end1=data$end1, start2=data$start2, end2 = data$end2, col = col_vec_comp[data$col],direction=1)
comparisonX = comparison(df)
return(comparisonX)
}

dna_segs = list(
make_dnaseq_eachID("hg38"),
make_dnaseq_eachID("aRhiBiv1"),
make_dnaseq_eachID("bCalAnn1"),
make_dnaseq_eachID("bStrHab1"),
make_dnaseq_eachID("bTaeGut2"),
make_dnaseq_eachID("fAnaTes1"),
make_dnaseq_eachID("fArcCen1"),
make_dnaseq_eachID("fAstCal1"),
make_dnaseq_eachID("fCotGob3"),
make_dnaseq_eachID("fGouWil2"),
make_dnaseq_eachID("fMasArm1"),
make_dnaseq_eachID("mLynCan4"),
make_dnaseq_eachID("mOrnAna1"),
make_dnaseq_eachID("mPhyDis1"),
make_dnaseq_eachID("mRhiFer1"),
make_dnaseq_eachID("rGopEvg1"),
make_dnaseq_eachID("sAmbRad1")
)
names(dna_segs)=sID_list

comparisons = list(
make_comp_eachID("01"),
make_comp_eachID("02"),
make_comp_eachID("03"),
make_comp_eachID("04"),
make_comp_eachID("05"),
make_comp_eachID("06"),
make_comp_eachID("07"),
make_comp_eachID("08"),
make_comp_eachID("09"),
make_comp_eachID("10"),
make_comp_eachID("11"),
make_comp_eachID("12"),
make_comp_eachID("13"),
make_comp_eachID("14"),
make_comp_eachID("15"),
make_comp_eachID("16")
)

pdf(paste0(wPATH,"output/BUSCO_genoPlotR_output/BUSCO_genoPlotR_colorred_by_",RefID,".pdf"),6,9)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,main=paste0("CoLink_Colorred_by_",RefID))
dev.off()

RefID = "fGouWil2"
iCNT_RefChr = 33
wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"
col_vec_comp = rainbow(iCNT_RefChr, alpha = 0.2)
col_vec_dnaseg = rainbow(iCNT_RefChr, alpha = 1) 

library(genoPlotR)
setwd(wPATH)

iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID



                                           
# function
# generate_dnaseg
make_dnaseq_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/dnaseg/dnaseg_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
if(sID==RefID){
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=col_vec_dnaseg[data$fill])
}else{
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=data$col, fill=data$fill)
}
dna_segX <- dna_seg(df,gene_type = "blocks")
return(dna_segX)
}

# generate_comp
make_comp_eachID <- function(compNumb){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/comparison/comparison",compNumb,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(start1=data$start1, end1=data$end1, start2=data$start2, end2 = data$end2, col = col_vec_comp[data$col],direction=1)
comparisonX = comparison(df)
return(comparisonX)
}

dna_segs = list(
make_dnaseq_eachID("hg38"),
make_dnaseq_eachID("aRhiBiv1"),
make_dnaseq_eachID("bCalAnn1"),
make_dnaseq_eachID("bStrHab1"),
make_dnaseq_eachID("bTaeGut2"),
make_dnaseq_eachID("fAnaTes1"),
make_dnaseq_eachID("fArcCen1"),
make_dnaseq_eachID("fAstCal1"),
make_dnaseq_eachID("fCotGob3"),
make_dnaseq_eachID("fGouWil2"),
make_dnaseq_eachID("fMasArm1"),
make_dnaseq_eachID("mLynCan4"),
make_dnaseq_eachID("mOrnAna1"),
make_dnaseq_eachID("mPhyDis1"),
make_dnaseq_eachID("mRhiFer1"),
make_dnaseq_eachID("rGopEvg1"),
make_dnaseq_eachID("sAmbRad1")
)
names(dna_segs)=sID_list

comparisons = list(
make_comp_eachID("01"),
make_comp_eachID("02"),
make_comp_eachID("03"),
make_comp_eachID("04"),
make_comp_eachID("05"),
make_comp_eachID("06"),
make_comp_eachID("07"),
make_comp_eachID("08"),
make_comp_eachID("09"),
make_comp_eachID("10"),
make_comp_eachID("11"),
make_comp_eachID("12"),
make_comp_eachID("13"),
make_comp_eachID("14"),
make_comp_eachID("15"),
make_comp_eachID("16")
)

pdf(paste0(wPATH,"output/BUSCO_genoPlotR_output/BUSCO_genoPlotR_colorred_by_",RefID,".pdf"),6,9)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,main=paste0("CoLink_Colorred_by_",RefID))
dev.off()

