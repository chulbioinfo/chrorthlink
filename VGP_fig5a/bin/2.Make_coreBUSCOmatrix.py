# This script to make input for both of BUSCO_circos plot and BUSCO_genoPlotR plot
# Version: 20200909
# Written by Chul Lee, SNU, Republic of Korea
# E-mail: chul.bioinfo@gmail.com

########################
#    Load_Libraries    #
########################
import sys, os, glob, math
import pandas as pd

#########################
#    GLOBAL_VARIABLE    #
#########################

Ref_BUSCO = "BUSCO_hg38.txt" # BUSCO_[assemblyID].txt
wdPATH = "../work/"
BUSCO_iPATH = wdPATH+"input/BUSCO/"

# HARDCODING to set order of sID lists
#sID_LIST = ["mLynCan4","mRhiFer1","mPhyDis1","mOrnAna1","bTaeGut1","bTaeGut2","bTaeGut2mat","bTaeGut2pat","bCalAnn1","bStrHab1","rGopEvg1","aRhiBiv1","sAmbRad1","fGouWil2","fArcCen1","fAstCal1","fMasArm1","fAnaTes1","fCotGob3"]

##########################
#    Default_VARIABLE    #
##########################
try:
  os.system("mkdir "+wdPATH )
except:
  pass
# file name : BUSCO_[assemblyID].txt
# Example header
'''
# BUSCO version is: 3.0.2 
# The lineage dataset is: vertebrata_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 2586)
# To reproduce this run: python /software/vertres/bin-external/busco-3.02/scripts/run_BUSCO.py -i aRhiBiv1.pri.cur.20190409.fasta -o aRhiBiv1.pri.cur.20190409.fasta.vert.human.again.busco -l /lustre/scratch116/vr/projects/vgp/user/mu2/busco_VGP/vertebrata_odb9/ -m genome -c 30 -t aRhiBiv1.pri.cur.20190409.fasta.vert.human.again.tmp/ -sp human
#
# Busco id	Status	Contig	Start	End	Score	Length
EOG090B000N	Missing
EOG090B0021	Fragmented	Super_scaffold_5	100847740	100879803	151.5	99
EOG090B002I	Fragmented	Super_scaffold_4	582304528	582383422	386.8	475
EOG090B0039	Missing
EOG090B003S	Fragmented	Super_scaffold_7	186107880	186122454	669.5	347
EOG090B004D	Fragmented	Super_scaffold_7	161565064	161600713	1060.7	535
EOG090B004X	Complete	Super_scaffold_11	107681758	107902650	4610.7	2270
'''
chrsize_iPATH = wdPATH+"input/chrsize/"
# file name : chrsize_[assemblyID].txt
# Example header
'''
Super_scaffold_1	595904407
Super_scaffold_10	103407986
Super_scaffold_11	128291458
Super_scaffold_12	83996221
Super_scaffold_13	80992113
Super_scaffold_14	82662308
Super_scaffold_15	74781573
'''
oNAME_full = wdPATH+"Full_Matrix_BUSCO.csv"
oNAME_core = wdPATH+"Core_Matrix_BUSCO.csv"


def Generate_BUSCO_matrix(BUSCO_iPATH):
  # Reading_data and Generating_matrix
  flist = glob.glob(BUSCO_iPATH + "BUSCO_*.txt")
  iCNT_Asm = len(flist)
  nAsmID_list = []
  nBuscoID_nInfo_dic = {}
  iCNT_order_Asm = 0 # start with one-based list
  for fNAME in flist:
    iCNT_order_Asm += 1
    AsmID = os.path.basename(fNAME).split(".")[0]
    nAsmID_list.append(AsmID)
    fpin = open(fNAME,'r')
    for line in fpin:
      if not line[0] == "#":
        part = line.strip("\n").split("\t")
        # Setting_default_BUSCO_info_list
        nBuscoID = part[0]
        if not nBuscoID in nBuscoID_nInfo_dic.keys():
          nBuscoID_nInfo_dic[nBuscoID] = []
          for i in range(iCNT_Asm):
            nBuscoID_nInfo_dic[nBuscoID].append([])
        # Generating_matrix
        nStatusInfo = part[1]
        if not nStatusInfo == "Missing":
          nScaffold = part[2]
          nStart = part[3]
          nEnd = part[4]
          nInfo = nStatusInfo+"&"+nScaffold+"&"+nStart+"&"+nEnd
          nBuscoID_nInfo_dic[nBuscoID][iCNT_order_Asm-1].append(nInfo)
    fpin.close()
  # Writing_matrix
  fpout_full = open(oNAME_full,'w')
  fpout_core = open(oNAME_core,'w')
  nHeader_full = "RefChrom"+','+"RefStartPos"+','+ "BuscoID"+','+','.join(nAsmID_list)+'\n'
  nHeader_core = "RefChrom"+','+"RefStartPos"+','+ "BuscoID"+','+','.join(nAsmID_list)+'\n'
  fpout_full.write(nHeader_full)
  fpout_core.write(nHeader_core)
  nRefAsmID = Ref_BUSCO.split(".")[0]
  nHeader_list = ["RefChrom","RefStartPos","BuscoID"]+nAsmID_list
  Index_RefAsmID  = nAsmID_list.index(nRefAsmID)
  for nBuscoID in nBuscoID_nInfo_dic.keys():
    tmpline = nBuscoID
    iCNT_SingleBUSCO = 0
    for nInfo_list in nBuscoID_nInfo_dic[nBuscoID]:
      if len(nInfo_list) > 0:
        tmpline += ','+'-'.join(nInfo_list)
        if nInfo_list[0].split("&")[0]=="Complete":
          iCNT_SingleBUSCO += 1
      else:
        tmpline += ','+''
    if len(nBuscoID_nInfo_dic[nBuscoID][Index_RefAsmID])>0:
      Ref_info = nBuscoID_nInfo_dic[nBuscoID][Index_RefAsmID][0]
      Ref_info_list = Ref_info.split("&")
      Ref_chr = Ref_info_list[1]
      Ref_start = Ref_info_list[2]
    else:
      Ref_chr = 'chr0'
      Ref_start = '0'
    tmpline = Ref_chr+','+Ref_start+','+tmpline+'\n'
    fpout_full.write(tmpline)
    if iCNT_SingleBUSCO == iCNT_Asm:
      fpout_core.write(tmpline)  
  fpout_full.close()  
  fpout_core.close()
Generate_BUSCO_matrix(BUSCO_iPATH)


# Sorting_parsing_complete_matrix
def ParseCompleteBUSCO_and_SortMatrix(fNAME):
  df = pd.read_csv(fNAME)
  df_sorted_by_StartPos = df.sort_values(by='RefStartPos' ,ascending=True)
  df_sorted_by_StartPos.to_csv(fNAME, header=True)
  fpin = open(fNAME,'r')
  line = fpin.readline()
  sID_list = line.strip('\n').split(",")[4:]
  nRefAsmID = Ref_BUSCO.split(".")[0]
  Index_RefAsmID  = sID_list.index(nRefAsmID)
  tmp_sID_list = [nRefAsmID]
  for i in range(len(sID_list)):
    if not i == Index_RefAsmID:
      tmp_sID_list.append(sID_list[i])
  nHeader = "BuscoID"+','+','.join(tmp_sID_list)+'\n'

  nChrom_nInfo_dic = {}
  for line in  fpin:
    part = line.strip('\n').split(",")[1:]
    if "chr" in part[0]:
      nChrom = part[0].split('chr')[1]
    else:
      nChrom = part[0]
    
    nStartPos = part[1]
    nBuscoID  = part[2]
    nBusco_info = part[3:]
    Ref_BUSCO_info = '&'.join(nBusco_info[Index_RefAsmID].split('&')[1:3])
    tmp_BuscoInfo_list = [Ref_BUSCO_info]
    for i in range(len(sID_list)):
      if not i == Index_RefAsmID:
        if not nBusco_info[i] == "":
          tmp_BUSCO_info_list = nBusco_info[i].split('&')
          nBuscoType = tmp_BUSCO_info_list[0]
          nBuscoChr  = tmp_BUSCO_info_list[1]
          nBuscoPos  = tmp_BUSCO_info_list[2]
          if nBuscoType == "Complete":
            tmp_BUSCO_info = nBuscoChr+'&'+nBuscoPos
            tmp_BuscoInfo_list.append(tmp_BUSCO_info)
          else:
            tmp_BuscoInfo_list.append("")
        else:
          tmp_BuscoInfo_list.append("")
      else:
        pass
    nInfo = nBuscoID+','+','.join(tmp_BuscoInfo_list)+'\n'
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
  fpout = open(fNAME,'w')
  fpout.write(nHeader)
  for nChrom in nChrom_list:
    for tmpline in nChrom_nInfo_dic[str(nChrom)]:
      fpout.write(tmpline)
  fpout.close()
ParseCompleteBUSCO_and_SortMatrix(oNAME_full)
ParseCompleteBUSCO_and_SortMatrix(oNAME_core)



