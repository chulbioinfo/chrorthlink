
library(genoPlotR)
setwd(wPATH)

iNAME_sID = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/sID_list.txt")
data = read.table(iNAME_sID,header=T)
sID_list = data$speciesID
sID = RefID



                                           
# function
# generate_dnaseg
make_dnaseq_eachID <- function(sID){
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/normalized/dnaseg/dnaseg_",sID,'.txt')
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
fNAME = paste0(wPATH,"BUSCO_genoPlotR_input/",RefID,"/normalized/comparison/comparison",compNumb,'.txt')
data = read.table(fNAME,sep="\t",header=T)
df = data.frame(start1=data$start1, end1=data$end1, start2=data$start2, end2 = data$end2, col = col_vec_comp[data$col],direction=1)
comparisonX = comparison(df)
return(comparisonX)
}

dna_segs = list(
