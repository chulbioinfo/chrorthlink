# This script to make input for both of BUSCO_circos plot and BUSCO_genoPlotR plot
# Version: 20210329
# Written by Chul Lee, SNU, Republic of Korea
# E-mail: chul.bioinfo@gmail.com

########################
#    Load libraries    #
########################
import sys
import os
import glob
import math
import pandas as pd


#########################
#    GLOBAL_VARIABLE    #
#########################
RefID_list = ["hg38"]#"hg38", "mOrnAna1", "bTaeGut2", "rGopEvg1", "aRhiBiv1", "sAmbRad1", "fGouWil2"]

sID_LIST = ["hg38","mLynCan4","mRhiFer1","mPhyDis1","mOrnAna1","rGopEvg1","bCalAnn1","bTaeGut1","bStrHab1","aRhiBiv1","fGouWil2","fAstCal1","fArcCen1","fCotGob3","fMasArm1","fAnaTes1","sAmbRad1"]

target_chr_name = "6" # or All

wPATH = "../work/"


##########################
#   Default variables    #
##########################

opt_sorting_chromosome_by_link = False # "chromosome_number" # or # "link"
opt_standardization_on = True # standardized genome sizes by the genome size of species with the max size of sum of chromosomes
opt_annotation_on = False # annoations on genoPlotR

figure_w = 5.4 # inch
figure_h = 6 # inch

fAlpha = 0.1
nRot_annot = "45"
iPATH_chrsize = wPATH+"input/chrsize/"




def reordering_matrix_by_sID_list(fNAME, sID_LIST):
  
  oNAME = "../work/reordered_"+os.path.basename(fNAME)
  fpout = open(oNAME,'w')
  nHeader = "BuscoID"
  for sID in sID_LIST:
    nHeader+=",BUSCO_"+sID
  fpout.write(nHeader+'\n')

  fpin = open(fNAME,'r')
  line = fpin.readline()
  part = line.strip().split(",")
  input_list = []
  for BUSCO_sID in part[1:]:
    sID = BUSCO_sID.split("_")[1]
    input_list.append(sID)

  for line in fpin:
    part = line.strip("\n").split(",")
    BUSCOname = part[0]
    info_list = part[1:]
    tmpline = BUSCOname
    for i in range(len(info_list)):
      tmpline += ','+info_list[input_list.index(sID_LIST[i])]
    fpout.write(tmpline+'\n')
  fpin.close()
  fpout.close()

reordering_matrix_by_sID_list("../work/Core_Matrix_BUSCO.csv",sID_LIST)
reordering_matrix_by_sID_list("../work/Full_Matrix_BUSCO.csv",sID_LIST)

#nColorPalette = "rainbow"
# better to re-order the columns of BUSCO results depend on the order you want to show
try:
  fNAME_BuscoCoreMatrix = wPATH+"reordered_Core_Matrix_BUSCO.csv"
  fpin = open(fNAME_BuscoCoreMatrix,'r')
  fpin.close()
except:
  try:
    fNAME_BuscoCoreMatrix = wPATH+"Core_Matrix_BUSCO.csv"
    fpin = open(fNAME_BuscoCoreMatrix,'r')
    fpin.close()
  except:
    print("Error - Core_Matrix_BUSCO.csv")
    sys.exit()


##################
#   Functions    #
##################
# Sorting_parsing_complete_matrix
def SortCompleteBuscoMatrix_By_RefID(RefID, fNAME):
  oNAME = os.path.basename(fNAME)
  oNAME = oNAME.replace(".csv","_sorted.csv")
  oNAME = wPATH+ oNAME

  fpout = open(oNAME,'w')
  fpin = open(fNAME,'r')
  nHeader = fpin.readline()
  part = nHeader.strip('\n').split(",")
  new_header = part[0]+','+'RefChr'+','+'RefStartPos'+','+','.join(part[1:])+'\n'
  fpout.write(new_header)
  
  sID_list = []
  for BUSCO_sID in part[1:]:
    sID = BUSCO_sID.split("BUSCO_")[1]
    sID_list.append(sID)
  
  for line in fpin:
    part= line.strip('\n').split(',')
    RefInfo = part[sID_list.index(RefID)+1]
    RefChrom = RefInfo.split("&")[0]
    RefStartPos = RefInfo.split("&")[1]
    newline = part[0]+','+RefChrom+','+RefStartPos+','+','.join(part[1:])+'\n'
    fpout.write(newline)
  fpout.close()
  
  df = pd.read_csv(oNAME)
  df_sorted_by_StartPos = df.sort_values(by='RefStartPos' ,ascending=True)
  df_sorted_by_StartPos.to_csv(oNAME, header=True)
  
  fpin = open(oNAME,'r')
  line = fpin.readline()
  
  nHeader = "BuscoID"
  for sID in sID_list:
    nHeader += ',' + "BUSCO_"+sID
  nHeader += '\n'

  nChrom_nInfo_dic = {}
  for line in fpin:
    part = line.strip('\n').split(",")[1:]
    nInfo = part[0]+','+','.join(part[3:])+'\n'
    nChrom = part[1]
    if "chr" in nChrom:
      nChrom = nChrom.split('chr')[1]
    else:
      nChrom = nChrom
    nChrom_nInfo_dic.setdefault(nChrom,[])
    nChrom_nInfo_dic[nChrom].append(nInfo)
  fpin.close()
  nChrom_int_list = []
  nChrom_str_list = []
  for nChrom in nChrom_nInfo_dic.keys():
    try:
      nChrom_int_list.append(int(nChrom))
    except:
      nChrom_str_list.append(str(nChrom))
  nChrom_int_list.sort()
  nChrom_str_list.sort()
  nChrom_list = nChrom_int_list + nChrom_str_list
  fpout = open(oNAME,'w')
  fpout.write(nHeader)
  for nChrom in nChrom_list:
    for tmpline in nChrom_nInfo_dic[str(nChrom)]:
      fpout.write(tmpline)
  fpout.close()
  return(oNAME)


def make_sID_SortedChrByBuscoOrder_dic(fNAME_BUSCO):
  fpin = open(fNAME_BUSCO,'r')
  nHeader = fpin.readline() # skip header
  fNAME_BUSCO_filetype = os.path.basename(fNAME_BUSCO).split(".")[-1] # .csv or .txt as tap
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
  iCNT_sID = len(sID_list)
 
  iCNT_nBuscoID = 0
  sID_ChrID_BuscoOrderList_dic = {}
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
      iCNT_nBuscoID += 1
      for i in range(1,len(part)):
        sID = sID_list[i-1]
        nChrom = part[i].split("&")[0]
        sID_ChrID_BuscoOrderList_dic.setdefault(sID,{})
        sID_ChrID_BuscoOrderList_dic[sID].setdefault(nChrom,[])
        sID_ChrID_BuscoOrderList_dic[sID][nChrom].append(iCNT_nBuscoID)
  fpin.close()

  sID_SortedChrByBuscoOrder_dic = {}
  for sID in sID_list:
    fMeanOrders_nChromList_dic = {}
    for nChrom in sID_ChrID_BuscoOrderList_dic[sID].keys():
      fMeanOrders = 0.0
      for iBuscoOrder in sID_ChrID_BuscoOrderList_dic[sID][nChrom]:
        fMeanOrders += iBuscoOrder
      fMeanOrders = fMeanOrders/len(sID_ChrID_BuscoOrderList_dic[sID][nChrom])
      fMeanOrders_nChromList_dic.setdefault(fMeanOrders,[])
      fMeanOrders_nChromList_dic[fMeanOrders].append(nChrom)

    fMeanOrders_list = []
    for fMeanOrders in fMeanOrders_nChromList_dic.keys():
      fMeanOrders_list.append(fMeanOrders)
    fMeanOrders_list.sort()
    
    SortedChrByBuscoOrder_list = []
    for fMeanOrders in fMeanOrders_list:
      if len(fMeanOrders_nChromList_dic[fMeanOrders])>1:
        MinBuscoOrder_list = []
        MinBuscoOrder_nChrom_dic = {}
        for nChrom in fMeanOrders_nChromList_dic[fMeanOrders]:
          MinBuscoOrder = sID_ChrID_BuscoOrderList_dic[sID][nChrom][0]
          MinBuscoOrder_list.append(MinBuscoOrder)
          MinBuscoOrder_nChrom_dic.setdefault(MinBuscoOrder,nChrom)
        MinBuscoOrder_list.sort()
        for MinBuscoOrder in MinBuscoOrder_list:
          nChrom = MinBuscoOrder_nChrom_dic[MinBuscoOrder]
          SortedChrByBuscoOrder_list.append(nChrom)
      elif len(fMeanOrders_nChromList_dic[fMeanOrders])==1:
        nChrom = fMeanOrders_nChromList_dic[fMeanOrders][0]
        SortedChrByBuscoOrder_list.append(nChrom)
    sID_SortedChrByBuscoOrder_dic.setdefault(sID,SortedChrByBuscoOrder_list)
  return(sID_SortedChrByBuscoOrder_dic)


def make_sID_SortedChrByChromNumb_dic(fNAME_BUSCO):
  fpin = open(fNAME_BUSCO,'r')
  nHeader = fpin.readline() # skip header
  fNAME_BUSCO_filetype = os.path.basename(fNAME_BUSCO).split(".")[-1] # .csv or .txt as tap
  if fNAME_BUSCO_filetype == "csv":
    part = nHeader.strip('\n').split(',')
  elif fNAME_BUSCO_filetype == "tsv":
    part = nHeader.strip('\n').split('\t')
  else:
    print("# BUSCO matrix file format error: ", fNAME_BUSCO)
    sys.exit()

  sID_nChrom_dic = {}  
  sID_list = []
  for BUSCO_sID in part[1:]:
    sID = BUSCO_sID.split("BUSCO_")[1]
    sID_list.append(sID)
    sID_nChrom_dic.setdefault(sID,{})
  iCNT_sID = len(sID_list)
  
  for line in fpin:
    part = line.strip('\n').split(",")
    for i in range(iCNT_sID):
      sID = sID_list[i]
      nChrom = part[i+1].split("&")[0]
      sID_nChrom_dic[sID].setdefault(nChrom,'')
  fpin.close()

  sID_SortedChrByChromNumb_dic = {}
  for sID in sID_list:
    nChrom_int_list = []
    nChrom_str_list = []
    nChrom_scf_list = []
    for nChrom in sID_nChrom_dic[sID].keys():
      if "chr" in nChrom:
        if not "_" in nChrom:
          nChrom = nChrom.split('chr')[1]
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
    for nChrom in nChrom_str_list:
      nChrom = "chr"+str(nChrom)
      nChrom_list.append(nChrom)
    sID_SortedChrByChromNumb_dic.setdefault(sID,nChrom_list)
  return(sID_SortedChrByChromNumb_dic)


def make_dnasegFiles_AdjLenDic_LinkColDic(RefID,iPATH,sID_SortedChrByBuscoOrder_dic,oPATH):
  maxLenID = ""
  maxLen   = 0
  sID_sumLen_dic = {}

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

    sumLen = 0.0
    for nChrom in nChrom_nLen_dic.keys():
      sumLen += int(nChrom_nLen_dic[nChrom])
    if sumLen > maxLen:
      maxLenID = sID
      maxLen = sumLen
    sID_sumLen_dic[sID]= sumLen

    
  sID_nChrom_iAddLen_dic = {}
  nRefChrom_nColCode_dic = {}
  
  flist = glob.glob(iPATH+"*.txt")
  for fNAME in flist:
    fpin =open(fNAME,'r')
    sID = os.path.basename(fNAME).split(".txt")[0].split("chrsize_")[1]
    SortedChr_list = sID_SortedChrByBuscoOrder_dic[sID]
    
    sumLen = sID_sumLen_dic[sID]
    if opt_standardization_on == True:
      normalize_value = maxLen/sumLen
    else:
      normalize_value = 1
    
    nChrom_nLen_dic = {}
    for line in fpin:
      part = line.strip('\n').split("\t")
      nChrom = part[0]
      nLen = part[1]
      if "chr"in nChrom:
        if not "_" in nChrom:
          nChrom_nLen_dic.setdefault(nChrom,nLen)
          nChrom_nLen_dic[nChrom]= str(int(int(nChrom_nLen_dic[nChrom])*normalize_value)) # nrmalization
    fpin.close()

    UnsortedChr_list = []
    for nChrom in nChrom_nLen_dic.keys():
      if not nChrom in SortedChr_list:
        UnsortedChr_list.append(nChrom)
    tmp_chr_list = SortedChr_list + UnsortedChr_list


    nChrom_int_list = []
    nChrom_str_list = []
    for nChrom in tmp_chr_list:
      nChrom = nChrom.split("chr")[1]
      try:
        nChrom_int_list.append(int(nChrom))
      except:
        nChrom_str_list.append(str(nChrom))

    nChrom_list = []
    for nChrom in tmp_chr_list:
      nChrom_list.append(nChrom)
    #if len(nChrom_str_list) == 0:
    #  nChrom_list.append("chrQ")

       
    sID_nChrom_iAddLen_dic.setdefault(sID,{})
    nChrom_iGenomicStart_iGenomicEnd_dic = {}
    iAddLen = 0
    iStart = 0
    iEnd = 0
    for nChrom in nChrom_list:
      sID_nChrom_iAddLen_dic[sID].setdefault(nChrom,iAddLen) # for comparison file
      iGenomicStart = iAddLen + 1
      iGenomicEnd = iAddLen + int(nChrom_nLen_dic[nChrom])
      nChrom_iGenomicStart_iGenomicEnd_dic.setdefault(nChrom,[iGenomicStart,iGenomicEnd])
      iAddLen = iGenomicEnd
      
    fpout_dnaseg = open(oPATH+"dnaseg/dnaseg_"+sID+'.txt','w')
    nHeader_dnaseg = "name"+'\t'+ "start" +'\t'+ "end" +'\t'+ "strand" +'\t'+ "col" +'\t'+ "fill" + '\t'+"lwd"+'\n'
    fpout_dnaseg.write(nHeader_dnaseg)
    
    fpout_annot = open(oPATH+"annotation/annotation_"+sID+'.txt','w')
    nHeader_annot = "x1"+'\t'+ "text" +'\t'+ "color" +'\t'+ "rot" +'\n'
    fpout_annot.write(nHeader_annot)
    if sID == RefID:
      iCNT = 5 # 1: black, 2: lightgrey, 3: grey, 4: white, 5: blue
      for nChrom in nChrom_list:
        iCNT +=1
        nName = nChrom.split("chr")[1]
        nStart = str(nChrom_iGenomicStart_iGenomicEnd_dic[nChrom][0])
        nEnd = str(nChrom_iGenomicStart_iGenomicEnd_dic[nChrom][1])
        nStrand = "1"
        nCol = "1"
        nFill = str(iCNT) # color code start with 5
        try:
          if nName == target_chr_name:
            nFill = "5"
        except:
          pass
        nLwd = "0.1"
        nRefChrom_nColCode_dic.setdefault(nChrom,nFill)
        tmpline_dnaseg = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ nCol +'\t'+ nFill +'\t'+ nLwd+ '\n'
        fpout_dnaseg.write(tmpline_dnaseg)
        tmpline_annot = str(int((int(nEnd)-int(nStart))/2.0)+int(nStart)) +'\t'+ nName +'\t'+ "black" +'\t'+ nRot_annot+ '\n'
        fpout_annot.write(tmpline_annot)
    else:
      iCNT = 5
      for nChrom in nChrom_list:
        iCNT +=1
        nName = nChrom.split("chr")[1]
        nStart = str(nChrom_iGenomicStart_iGenomicEnd_dic[nChrom][0])
        nEnd = str(nChrom_iGenomicStart_iGenomicEnd_dic[nChrom][1])
        nStrand = "1"
        nLwd = "0.1"
        nCol = "1"
        if not int(iCNT/2.0)*2==iCNT: # odd number
          nFill = "4"
        else:# even number
          nFill = "4"
        tmpline_dnaseg = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ nCol +'\t'+ nFill +'\t'+ nLwd+ '\n'
        fpout_dnaseg.write(tmpline_dnaseg)
        tmpline_annot = str(int((int(nEnd)-int(nStart))/2.0)+int(nStart)) +'\t'+ nName +'\t'+ "black" +'\t'+ nRot_annot+ '\n'
        fpout_annot.write(tmpline_annot)
    fpout_dnaseg.close()
    fpout_annot.close()
  return sID_nChrom_iAddLen_dic, nRefChrom_nColCode_dic, sID_sumLen_dic

def make_comparisonFiles(RefID, fNAME_BUSCO, sID_nChrom_iAddLen_dic, nRefChrom_nColCode_dic, sID_sumLen_dic, oPATH):
  LinkColCode_list = []

  maxLen = 0
  for sID in sID_sumLen_dic.keys():
    if sID_sumLen_dic[sID] > maxLen:
      maxLen = sID_sumLen_dic[sID]
  
  iCNT_sID = len(sID_nChrom_iAddLen_dic)
  sID_nPos_dic = {}
  fpin = open(fNAME_BUSCO,'r')
  nHeader = fpin.readline()
  
  part = nHeader.strip('\n').split(',')
  sID_list = []
  for BUSCO_sID in part[1:]:
    sID_list.append(BUSCO_sID.split("BUSCO_")[1])
  for line in fpin:
    part = line.strip('\n').split(',')
    iCNT_chr = 0
    for nInfo in part[1:]:
      if "chr" in nInfo:
        if not "_" in nInfo:
          iCNT_chr += 1
    if iCNT_chr == iCNT_sID:
      for i in range(1,len(part)):
        sID = sID_list[i-1]
        if opt_standardization_on == True:
          normalized_value = maxLen/sID_sumLen_dic[sID]
        else:
          normalized_value = 1
        nChrom = part[i].split("&")[0]
        iPos_Chr = int(int(part[i].split("&")[1])*normalized_value)
        iAdjLen = int(sID_nChrom_iAddLen_dic[sID][nChrom])
        iPos_Genome = iAdjLen + iPos_Chr
        sID_nPos_dic.setdefault(sID,[])
        sID_nPos_dic[sID].append(str(iPos_Genome))
        if sID == RefID:
          if  nChrom == "chr"+target_chr_name:
            LinkColCode_list.append("5")
          else:
            LinkColCode_list.append(nRefChrom_nColCode_dic[nChrom])
  fpin.close()

  # Write links on dnaseg
  for sID in sID_list:
    lines_fill = ''
    lines_back = ''
    lines_line = ''
    fpin = open(oPATH+"dnaseg/dnaseg_"+sID+'.txt','r')
    nHeader = fpin.readline()
    for line in fpin:
      part = line.strip('\n').split("\t") # tmpline = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ nCol +'\t'+ nFill +'\t'+ nLwd+ '\n'
      nName   = part[0]
      nStart  = part[1]
      nEnd    = part[2]
      nStrand = part[3]
      nCol    = part[4]
      nFill   = part[5]
      nLwd    = part[6]
      tmpline_back = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ "4" +'\t'+ "4" +'\t'+ nLwd+ '\n'
      lines_back += tmpline_back
      tmpline_fill = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ nCol +'\t'+ nFill +'\t'+ nLwd+ '\n'
      lines_fill += tmpline_fill
      tmpline_line = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ nCol +'\t'+ "" +'\t'+ "1"+ '\n'
      lines_line += tmpline_line
    fpin.close()
    
    fpout = open(oPATH+"dnaseg/dnaseg_"+sID+'.txt','w')
    fpout.write(nHeader)
    fpout.write(lines_back)
    fpout.write(lines_fill)
    # tmpline = nName+'\t'+ nStart +'\t'+ nEnd  +'\t'+ nStrand +'\t'+ nCol +'\t'+ nFill +'\t'+ nLwd+ '\n'
    for i in range(len(LinkColCode_list)):
      nName = "L"
      iPos = int(sID_nPos_dic[sID][i])
      nStrand = "1"
      nCol = LinkColCode_list[i]
      nFill = LinkColCode_list[i]
      nLwd = "0.01"
      tmpline = nName+'\t'+ str(iPos) +'\t'+ str(iPos) +'\t'+ nStrand +'\t'+ nCol +'\t'+ nFill +'\t'+ nLwd+ '\n'
      fpout.write(tmpline)
    fpout.write(lines_line)
    fpout.close()

  # Write comparison
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

    tmp_line_target = ''
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
 
########
# Main #
########

for RefID in RefID_list:
  oPATH = wPATH + "BUSCO_genoPlotR_input/"+RefID+"/"
  try:
    os.system("mkdir ./BUSCO_genoPlotR_input")
    os.system("mkdir "+oPATH)
    os.system("mkdir "+oPATH+"dnaseg")
    os.system("mkdir "+oPATH+"comaprison")
    os.system("mkdir "+oPATH+"annotation")
    os.system("mkdir ./output")

  except:
    pass

  # species name on PDF
  fpout = open(oPATH+"sID_list.txt",'w')
  fpout.write("speciesID\n")
  for sID in sID_LIST:
    fpout.write(sID+'\n')
  fpout.close()
    
  fNAME_SortedBuscoCoreMatrix = SortCompleteBuscoMatrix_By_RefID(RefID, fNAME_BuscoCoreMatrix)
  if opt_sorting_chromosome_by_link == True:
    sID_SortedChr_dic = make_sID_SortedChrByBuscoOrder_dic(fNAME_SortedBuscoCoreMatrix)
  else:
    sID_SortedChr_dic = make_sID_SortedChrByChromNumb_dic(fNAME_SortedBuscoCoreMatrix)
  sID_nChrom_iAddLen_dic, nRefChrom_nColCode_dic, sID_sumLen_dic = make_dnasegFiles_AdjLenDic_LinkColDic(RefID,iPATH_chrsize,sID_SortedChr_dic,oPATH)
  sID_list , iCNT_RefChr = make_comparisonFiles(RefID, fNAME_BuscoCoreMatrix, sID_nChrom_iAddLen_dic, nRefChrom_nColCode_dic, sID_sumLen_dic, oPATH)
  


