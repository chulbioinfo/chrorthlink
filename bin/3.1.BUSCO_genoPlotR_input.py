# This script to make input for both of BUSCO_circos plot and BUSCO_genoPlotR plot
# Version: 20200909
# Written by Chul Lee, SNU, Republic of Korea
# E-mail: chul.bioinfo@gmail.com

import sys, os, glob, math
import pandas as pd
#########################
#    GLOBAL_VARIABLE    #
#########################

# BUSCO_[assemblyID].txt
#RefID = "hg38" 
#RefID = "mOrnAna1" 
#RefID = "bTaeGut1" 
#RefID = "rGopEvg1" 
#RefID = "aRhiBiv1" 
#RefID = "sAmbRad1" 
#RefID = "fCotGob3"
fAlpha = 0.2

RefID_list = ["hg38", "mOrnAna1", "bTaeGut2", "rGopEvg1", "aRhiBiv1", "sAmbRad1", "fGouWil2"]

iPATH_chrsize = "../example/VGP_1st/work/input/chrsize/"
wPATH = "../example/VGP_1st/work/"

#fNAME_BuscoCoreMatrix = wPATH+"reordered_Core_Matrix_BUSCO.csv"# it is better to re-order the columns of BUSCO results depend on the order you want to show
fNAME_BuscoCoreMatrix = wPATH+"Core_Matrix_BUSCO.csv"

##################
#   Functions    #
##################
def make_dnasegFiles_AdjLenDic_LinkColDic(RefID,iPATH,oPATH):
  sID_nChrom_iAddLen_dic = {}
  nRefChrom_nColCode_dic = {}
  
  flist = glob.glob(iPATH+"*.txt")
  for fNAME in flist:
    fpin =open(fNAME,'r')
    sID = os.path.basename(fNAME).split(".txt")[0].split("chrsize_")[1]
    nChrom_nLen_dic = {}
    for line in fpin:
      part = line.strip('\n').split("\t")
      nChrom = part[0]
      nLen = part[1]
      if "chr"in nChrom:
        if not "_" in nChrom:
          nChrom_nLen_dic.setdefault(nChrom,nLen)
    fpin.close()
    
    nChrom_int_list = []
    nChrom_str_list = []
    for nChrom in nChrom_nLen_dic.keys():
      nChrom = nChrom.split("chr")[1]
      try:
        nChrom_int_list.append(int(nChrom))
      except:
        nChrom_str_list.append(str(nChrom))
    nChrom_int_list.sort()
    nChrom_str_list.sort()
    nChrom_list = []
    for nChrom in nChrom_int_list:
      nChrom = "chr"+str(nChrom)
      nChrom_list.append(nChrom)
    if len(nChrom_str_list) > 0:
      for nChrom in nChrom_str_list:
        nChrom = "chr"+str(nChrom)
        nChrom_list.append(nChrom)
    else:
      nChrom_list.append("chrQ")

    sID_nChrom_iAddLen_dic.setdefault(sID,{})
    nChrom_iGenomicStart_iGenomicEnd_dic = {}
    iAddLen = 0
    iStart = 0
    iEnd = 0
    for nChrom in nChrom_list:
      sID_nChrom_iAddLen_dic[sID].setdefault(nChrom,iAddLen) # for comparison file
      if not nChrom == "chrQ":
        iGenomicStart = iAddLen + 1
        iGenomicEnd = iAddLen + int(nChrom_nLen_dic[nChrom])
      else:
        iGenomicStart = iAddLen
        iGenomicEnd = iAddLen
      nChrom_iGenomicStart_iGenomicEnd_dic.setdefault(nChrom,[iGenomicStart,iGenomicEnd])
      iAddLen = iGenomicEnd

    fpout = open(oPATH+"dnaseg/dnaseg_"+sID+'.txt','w')
    nHeader = "name"+'\t'+ "start" +'\t'+ "end" +'\t'+ "strand" +'\t'+ "col" +'\t'+ "fill" + '\n'
    fpout.write(nHeader)
    if sID == RefID:
      iCNT = 0
      for nChrom in nChrom_list:
        iCNT +=1
        nName = nChrom.split("chr")[1]
        nStart = str(nChrom_iGenomicStart_iGenomicEnd_dic[nChrom][0])
        nEnd = str(nChrom_iGenomicStart_iGenomicEnd_dic[nChrom][1])
        nStrand = "1"
        nCol = "black"
        nFill = str(iCNT)
        nRefChrom_nColCode_dic.setdefault(nChrom,nFill)
        tmpline = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ nCol +'\t'+ nFill + '\n'
        fpout.write(tmpline)
    else:
      iCNT = 0
      for nChrom in nChrom_list:
        iCNT +=1
        nName = nChrom.split("chr")[1]
        nStart = str(nChrom_iGenomicStart_iGenomicEnd_dic[nChrom][0])
        nEnd = str(nChrom_iGenomicStart_iGenomicEnd_dic[nChrom][1])
        nStrand = "1"
        nCol = "black"
        if not int(iCNT/2.0)*2==iCNT: # odd number
          nFill = "grey"
        else:# even number
          nFill = "lightgrey"
        tmpline = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ nCol +'\t'+ nFill + '\n'
        fpout.write(tmpline)
    fpout.close()
  return sID_nChrom_iAddLen_dic, nRefChrom_nColCode_dic

def make_comparisonFiles(RefID, fNAME_BUSCO, sID_nChrom_iAddLen_dic, nRefChrom_nColCode_dic, oPATH):
  LinkColCode_list = []
  iCNT_sID = len(sID_nChrom_iAddLen_dic)
  sID_nPos_dic = {}
  fpin = open(fNAME_BUSCO,'r')
  fNAME_BUSCO_filetype = os.path.basename(fNAME_BUSCO).split(".")[-1] # .csv or .txt as tap
  nHeader = fpin.readline()
  
  if fNAME_BUSCO_filetype == "csv":
    part = nHeader.strip('\n').split(',')
  elif fNAME_BUSCO_filetype == "tsv":
    part = nHeader.strip('\n').split('\t')
  else:
    print("# BUSCO matrix file format error: ", fNAME_BUSCO)
    sys.exit()
    
  sID_list = []
  for BUSCO_sID in part[1:]:
    sID_list.append(BUSCO_sID.split("BUSCO_")[1])

  for line in fpin:
    if fNAME_BUSCO_filetype == "csv":
      part = line.strip('\n').split(',')
    elif fNAME_BUSCO_filetype == "txt":
      part = line.strip('\n').split('\t')
    iCNT_chr = 0
    for nInfo in part[1:]:
      if "chr" in nInfo:
        if not "_" in nInfo:
          iCNT_chr += 1
    if iCNT_chr == iCNT_sID:
      for i in range(1,len(part)):
        sID = sID_list[i-1]
        nChrom = part[i].split("&")[0]
        iPos_Chr = int(part[i].split("&")[1])
        iAdjLen = int(sID_nChrom_iAddLen_dic[sID][nChrom])
        iPos_Genome = iAdjLen + iPos_Chr
        sID_nPos_dic.setdefault(sID,[])
        sID_nPos_dic[sID].append(str(iPos_Genome))
        if sID == RefID:
          LinkColCode_list.append(nRefChrom_nColCode_dic[nChrom])
  fpin.close()

  iCNT_comparison = 0
  tmp_s1_Pos_list = sID_nPos_dic[sID_list[0]]
  for sID in sID_list[1:]:
    iCNT_comparison += 1
    if len(str(iCNT_comparison))==1:
      oNAME = oPATH+ "comparison/comparison0"+str(iCNT_comparison)+".csv"
    else:
      oNAME = oPATH+ "comparison/comparison"+str(iCNT_comparison)+".csv"
    fpout=open(oNAME,'w')
    nHeader = "start1,end1,start2,end2,col\n"
    fpout.write(nHeader)
    for i in range(len(LinkColCode_list)):
      s1_nPos = tmp_s1_Pos_list[i]
      s2_nPos = sID_nPos_dic[sID][i]
      tmpline = s1_nPos+','+s1_nPos+','+s2_nPos+','+s2_nPos+','+LinkColCode_list[i]+'\n'
      fpout.write(tmpline)
    fpout.close()
    tmp_s1_Pos_list = sID_nPos_dic[sID]

  # Reordering_by_LinkCol
  compNAME_list = glob.glob(oPATH+ "comparison/comparison*.csv")
  for compNAME in compNAME_list:
    
    df = pd.read_csv(compNAME)
    df_sorted_by_StartPos = df.sort_values(by='col',ascending=False)
    df_sorted_by_StartPos.to_csv(compNAME, header=True)
    fpin = open(compNAME,'r')
    oNAME = oPATH+ "comparison/"+os.path.basename(compNAME).split(".")[0]+".txt"
    fpout = open(oNAME,'w')
    for line in fpin:
      tmpline = '\t'.join(line.split(",")[1:])
      fpout.write(tmpline)
    fpin.close()
    fpout.close()


    
  fpout = open(oPATH+"sID_list.txt",'w')
  fpout.write("speciesID\n")
  for sID in sID_list:
    fpout.write(sID+'\n')
  fpout.close()
  iCNT_RefChr = len(nRefChrom_nColCode_dic.keys())
  return sID_list, iCNT_RefChr
 
# write_Rscript
def write_Rscript(RefID, sID_list , iCNT_RefChr):
 
  fpout = open(wPATH+"BUSCO_genoPlotR_input/3.1.BUSCO_genoPlotR_by_"+RefID+".R",'w')

  tmpline = "RefID = "+'"'+RefID+'"'+"\n"
  fpout.write(tmpline)

  tmpline = "iCNT_RefChr = "+str(iCNT_RefChr+int(iCNT_RefChr*0.4))+"\n"
  fpout.write(tmpline)

  #tmpline = "wPATH = "+'"'+wPATH+'"'+'\n' 
  tmpline = "wPATH = "+'"'+"E:/GoogleDrive/Research/2020/ChrOrthoLink/example/VGP_1st/work/"+'"'+'\n'
  fpout.write(tmpline)

  tmpline = "col_vec_comp = rainbow(iCNT_RefChr, alpha = "+ str(fAlpha) +")\n"
  fpout.write(tmpline)

  fpin = open("Rscript/tmplete_BUSCO_genoPlotR.R",'r')
  for line in fpin:
    fpout.write(line)
  fpin.close()

  for sID in sID_list[:-1]:
    tmpline = "make_dnaseq_eachID("+'"'+sID+'"'+"),\n"
    fpout.write(tmpline)
  tmpline = "make_dnaseq_eachID("+'"'+sID_list[-1]+'"'+")\n)\nnames(dna_segs)=sID_list\n\ncomparisons = list(\n"
  fpout.write(tmpline)
  for i in range(1,len(sID_list)-1):
    if len(str(i))==1:
      compID = "0"+str(i)
    else:
      compID = str(i)
    tmpline =  "make_comp_eachID("+'"'+compID+'"'+ "),\n"
    fpout.write(tmpline)
  if len(str(len(sID_list)-1))==1:
    compID = "0"+str(len(sID_list)-1)
  else:
    compID = str(len(sID_list)-1)
  tmpline =  "make_comp_eachID("+'"'+compID+'"'+ ")\n)\n\n"
  fpout.write(tmpline)

  tmpline = 'pdf(paste0(wPATH,'+'"'+'output/BUSCO_genoPlotR_output/BUSCO_genoPlotR_colorred_by_'+'"'+',RefID,'+'"'+'.pdf'+'"'+'),6,9)\n'
  fpout.write(tmpline)

  tmpline = 'plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,main=paste0('+'"'+'CoLink_Colorred_by_'+'"'+',RefID))\ndev.off()\n\n'
  fpout.write(tmpline)
  fpout.close()


########
# Main #
########

for RefID in RefID_list:
  oPATH = wPATH + "BUSCO_genoPlotR_input/"+RefID+"/"
  try:
    os.system("mkdir "+wPATH +"BUSCO_genoPlotR_input")
    os.system("mkdir "+oPATH)
    os.system("mkdir "+oPATH+"dnaseg")
    os.system("mkdir "+oPATH+"comparison")
  except:
    pass
  sID_nChrom_iAddLen_dic, nRefChrom_nColCode_dic = make_dnasegFiles_AdjLenDic_LinkColDic(RefID,iPATH_chrsize,oPATH)
  sID_list , iCNT_RefChr = make_comparisonFiles(RefID, fNAME_BuscoCoreMatrix, sID_nChrom_iAddLen_dic, nRefChrom_nColCode_dic, oPATH)
  write_Rscript(RefID,sID_list,iCNT_RefChr)

fpout = open("./4.1.BUSCO_genoPlotR_by_EachRefID.R",'w')
for RefID in RefID_list:
  fpin = open(wPATH+"BUSCO_genoPlotR_input/3.1.BUSCO_genoPlotR_by_"+RefID+".R",'r')
  for line in fpin:
    fpout.write(line)
  fpin.close()
fpout.close()
  
