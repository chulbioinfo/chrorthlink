# install.packages("RColorBrewer")
# install.packages("genoPlotR") 

##################
# Setting Colors #
##################

library(RColorBrewer)
makeTransparent = function(..., alpha=0.5) {
if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
alpha = floor(255*alpha)
newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
.makeTransparent = function(col, alpha) {
  rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
return(newColor)
}
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
iCNT_RefChr = 35
col_vec_dnaseg = c("black","grey","lightgrey","white","blue",makeTransparent(col_vector, alpha = 0.1), makeTransparent("blue", alpha = 0.1))
col_vec_comp = c("black","grey","lightgrey","white","blue",makeTransparent(col_vector, alpha = 0.1), makeTransparent("blue", alpha = 0.1))



###############################################################
# Visualization of chromosomal orthologous link : ChrOrthLink #
###############################################################
library(genoPlotR)

RefID = "hg38"
wPATH = "../work/"
#wPATH = "E:/GoogleDrive/Research/2020/ChrOrthoLink/VGP_fig5a/work/"
setwd(wPATH)

# read species IDs
iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID

                                           
# function
# generate_dnaseg
make_dnaseg_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/dnaseg/dnaseg_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(name=data$name, start=data$start, end=data$end, strand=data$strand, col=col_vec_dnaseg[data$col], fill=col_vec_dnaseg[data$fill], lwd = data$lwd)
dna_segX <- dna_seg(df,gene_type = "blocks")
return(dna_segX)
}

# generate_annot
make_annotation_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/annotation/annotation_",sID,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(x1=data$x1, text=data$text, color=data$color, rot = data$rot)
annotationX = as.annotation(df, x2 = NA)
return(annotationX)
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
make_dnaseg_eachID("hg38"),
make_dnaseg_eachID("mLynCan4"),
make_dnaseg_eachID("mRhiFer1"),
make_dnaseg_eachID("mPhyDis1"),
make_dnaseg_eachID("mOrnAna1"),
make_dnaseg_eachID("rGopEvg1"),
make_dnaseg_eachID("bCalAnn1"),
make_dnaseg_eachID("bTaeGut1"),
make_dnaseg_eachID("bStrHab1"),
make_dnaseg_eachID("aRhiBiv1"),
make_dnaseg_eachID("fGouWil2"),
make_dnaseg_eachID("fAstCal1"),
make_dnaseg_eachID("fArcCen1"),
make_dnaseg_eachID("fCotGob3"),
make_dnaseg_eachID("fMasArm1"),
make_dnaseg_eachID("fAnaTes1"),
make_dnaseg_eachID("sAmbRad1")
)
names(dna_segs)=sID_list

annotations = list(
make_annotation_eachID("hg38"),
make_annotation_eachID("mLynCan4"),
make_annotation_eachID("mRhiFer1"),
make_annotation_eachID("mPhyDis1"),
make_annotation_eachID("mOrnAna1"),
make_annotation_eachID("rGopEvg1"),
make_annotation_eachID("bCalAnn1"),
make_annotation_eachID("bTaeGut1"),
make_annotation_eachID("bStrHab1"),
make_annotation_eachID("aRhiBiv1"),
make_annotation_eachID("fGouWil2"),
make_annotation_eachID("fAstCal1"),
make_annotation_eachID("fArcCen1"),
make_annotation_eachID("fCotGob3"),
make_annotation_eachID("fMasArm1"),
make_annotation_eachID("fAnaTes1"),
make_annotation_eachID("sAmbRad1")
)
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

left_offset = rep(c(0), times = length(dna_segs))

pdf(paste0(wPATH,"output/VGP_Fig5a_(5.4x6).pdf"),5.4,6)
 plot_gene_map(dna_segs=dna_segs, comparisons=comparisons, minimum_gap_size = 0.5, offsets = left_offset, scale_cex = 0, scale = FALSE)
dev.off()

