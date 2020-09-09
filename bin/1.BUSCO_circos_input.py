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
wdPATH = "../example/VGP_1st/work/"
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
BuscoCircosInput_oPATH = wdPATH+"BUSCO_circos_input/"
try:
  os.system("mkdir "+BuscoCircosInput_iPATH)
except:
  pass
BuscoCircosInput_byChr_oPATH = wdPATH+"BUSCO_circos_input_byChr/"
try:
  os.system("mkdir "+BuscoCircosInput_byChr_oPATH)
except:
  pass

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


# Generate_LinkFile_from_BUSCO_matrix
def Generate_LinkFile_from_FullMatrix(oNAME_full):
  # Reading sorted matrix to generate_link_matrixs
  fpin = open(oNAME_full,'r')
  line = fpin.readline()
  nHeader_list = line.strip('\n').split(',')
  sID_link_dic = {}
  for line in fpin:
    part = line.strip('\n').split(',')
    for i in range(2,len(nHeader_list)):
      sID = nHeader_list[i].split("_")[1]
      if not part[1]=="":
        if not part[i]=="":
          Ref_info = '\t'.join(part[1].split("&"))+'\t'+part[0]
          Que_info = '\t'.join(part[i].split("&"))+'\t'+part[0]
          link_info = Ref_info + '\t' + Que_info +'\n'
          sID_link_dic.setdefault(sID,'')
          sID_link_dic[sID] += link_info
  fpin.close()
  # Generating_link_file
  for sID in sID_link_dic.keys():
    fpout = open(BuscoCircosInput_oPATH+"FullLink_"+sID+'.tsv','w')
    link_header = 'chr1\tpo1\tgene1\tchr2\tpo2\tgene2\n'
    fpout.write(link_header+sID_link_dic[sID])
    fpout.close()
  # Generating_link_chr
  RefChr_link_dic = {}
  link_flist = glob.glob(BuscoCircosInput_oPATH+"FullLink_*.tsv")
  for fNAME_link in link_flist:
    # reading link files and setting list of chr of each sID
    # Writing_link color by chr
    sID = os.path.basename(fNAME_link).split(".tsv")[0].split("Link_")[1]
    fpin = open(fNAME_link,'r')
    nHeader_LinkByChr = fpin.readline()
    for line in fpin:
      part = line.strip('\n').split("\t")
      Ref_chr = part[0]
      Ref_pos = part[1]
      BUSCOid = part[2]
      Que_chr = sID+"_"+part[3]
      Que_pos = part[4]
      tmpline_link = Ref_chr +'\t'+ Ref_pos +'\t'+ BUSCOid +'\t'+ Que_chr +'\t'+ Que_pos +'\t'+ sID+"_"+BUSCOid +'\n'
      RefChr_link_dic.setdefault(Ref_chr,'')
      RefChr_link_dic[Ref_chr]+=tmpline_link
    fpin.close()
  for RefChr in RefChr_link_dic.keys():
    fpout = open(BuscoCircosInput_byChr_oPATH+"FullLink_"+RefChr+'.tsv','w')
    fpout.write(nHeader_LinkByChr+RefChr_link_dic[RefChr])
    fpout.close()
Generate_LinkFile_from_FullMatrix(oNAME_full)


def Generate_LinkFile_from_CoreMatrix(oNAME_core):
  # Reading sorted matrix to generate_link_matrixs
  fpin = open(oNAME_core,'r')
  line = fpin.readline()
  nHeader_list = line.strip('\n').split(',')
  sID_link_dic = {}
  for line in fpin:
    part = line.strip('\n').split(',')
    for i in range(2,len(nHeader_list)):
      sID = nHeader_list[i].split("_")[1]
      if not part[1]=="":
        if not part[i]=="":
          Ref_info = '\t'.join(part[1].split("&"))+'\t'+part[0]
          Que_info = '\t'.join(part[i].split("&"))+'\t'+part[0]
          link_info = Ref_info + '\t' + Que_info +'\n'
          sID_link_dic.setdefault(sID,'')
          sID_link_dic[sID] += link_info
  fpin.close()
  # Generating_link_file
  for sID in sID_link_dic.keys():
    fpout = open(BuscoCircosInput_oPATH+"CoreLink_"+sID+'.tsv','w')
    link_header = 'chr1\tpo1\tgene1\tchr2\tpo2\tgene2\n'
    fpout.write(link_header+sID_link_dic[sID])
    fpout.close()
  # Generating_link_chr
  RefChr_link_dic = {}
  link_flist = glob.glob(BuscoCircosInput_oPATH+"CoreLink_*.tsv")
  for fNAME_link in link_flist:
    # reading link files and setting list of chr of each sID
    # Writing_link color by chr
    sID = os.path.basename(fNAME_link).split(".tsv")[0].split("Link_")[1]
    fpin = open(fNAME_link,'r')
    nHeader_LinkByChr = fpin.readline()
    for line in fpin:
      part = line.strip('\n').split("\t")
      Ref_chr = part[0]
      Ref_pos = part[1]
      BUSCOid = part[2]
      Que_chr = sID+"_"+part[3]
      Que_pos = part[4]
      tmpline_link = Ref_chr +'\t'+ Ref_pos +'\t'+ BUSCOid +'\t'+ Que_chr +'\t'+ Que_pos +'\t'+ sID+"_"+BUSCOid +'\n'
      RefChr_link_dic.setdefault(Ref_chr,'')
      RefChr_link_dic[Ref_chr]+=tmpline_link
    fpin.close()
  for RefChr in RefChr_link_dic.keys():
    fpout = open(BuscoCircosInput_byChr_oPATH+"CoreLink_"+RefChr+'.tsv','w')
    fpout.write(nHeader_LinkByChr+RefChr_link_dic[RefChr])
    fpout.close()
Generate_LinkFile_from_CoreMatrix(oNAME_core)


# Generating_chrsize_chrname_chrcolor_linkcolor_from_link
def Generating_chrsize_chrname_chrcolor_linkcolor_from_link(nLinkType):# ex) nLinkType = Full
  link_flist = glob.glob(BuscoCircosInput_byChr_oPATH+nLinkType+"Link_*.tsv")
  for fNAME_link in link_flist:
    # reading link files and setting list of chr of each sID
    # Writing_link color by chr
    ChrID = os.path.basename(fNAME_link).split(".tsv")[0].split(nLinkType+"Link_")[1]

    LinkColor_fNAME = BuscoCircosInput_byChr_oPATH+nLinkType+"LinkColor_"+ChrID+'.tsv'
    
    fpout_LinkColor = open(LinkColor_fNAME,'w')
    LinkColor_nHeader = "link_ID\tcolor\n"
    fpout_LinkColor.write(LinkColor_nHeader)

    fpin = open(fNAME_link,'r')
    fpin.readline()
    sID_list = []
    
    Ref_Chr_list = []
    Que_Chr_list = []
    iCNT_linkID = 0
    iColorCode_Link = 0
    for line in fpin:
      iCNT_linkID += 1
      nCNT_linkID = str(iCNT_linkID)
      part = line.strip('\n').split("\t")
      Que_chr = part[3]
      sID = part[3].split("_")[0]
      Ref_chr = part[0]
      if not Ref_chr in Ref_Chr_list:
        Ref_Chr_list.append(Ref_chr)
      if not sID in sID_list:
        sID_list.append(sID)
        iColorCode_Link += 1
        nColorCode_Link = str(iColorCode_Link)
      tmpline = nCNT_linkID + '\t' + nColorCode_Link + '\n'
      fpout_LinkColor.write(tmpline)
      if not Que_chr in Que_Chr_list:
        Que_Chr_list.append(Que_chr)
    fpin.close()
    fpout_LinkColor.close()
    

    # Reordering Que_Chr_list
    nChrom_int_dic = {}
    nChrom_str_dic = {}
    nScaffold_str_dic = {}
    for sID_nChrom in Que_Chr_list:
      sID = sID_nChrom.split("_")[0]
      nChrom_int_dic.setdefault(sID,[])
      nChrom_str_dic.setdefault(sID,[])
      nScaffold_str_dic.setdefault(sID,[])
      if "chr" in sID_nChrom:
        nChrom = sID_nChrom.split("chr")[1]
        try:
          nChrom_int_dic[sID].append(int(nChrom))
        except:
          nChrom_str_dic[sID].append(str(nChrom))
      else:
        nScaffold_str_dic[sID].append(sID_nChrom)
    tmp_sID_list = []
    for sID in nChrom_int_dic.keys():
      tmp_sID_list.append(sID)
      nChrom_int_dic[sID].sort()
      nChrom_str_dic[sID].sort()
      nScaffold_str_dic[sID].sort()
    tmp_sID_list.sort()
    
    Que_Chr_list = []
    for sID in tmp_sID_list:
      for nChrom in nChrom_int_dic[sID]:
        Que_Chr_list.append(sID+"_chr"+str(nChrom))
      for nChrom in nChrom_str_dic[sID]:
        Que_Chr_list.append(sID+"_chr"+str(nChrom))
      Que_Chr_list += nScaffold_str_dic[sID]


    # Reading_chrsize_Ref
    Ref_id = Ref_BUSCO.split("_")[1].split(".")[0]
    Ref_chrsize_fNAME = chrsize_iPATH+'chrsize_'+Ref_id+'.txt' # hard coding for chrsize files - need to be updated
    fpin = open(Ref_chrsize_fNAME,'r') 
    Ref_nChrom_iLen_dic = {}
    for line in fpin:
      part = line.strip().split("\t")
      Ref_nChrom_iLen_dic.setdefault(part[0],part[1])
    fpin.close()

    # Reading_chrsize_Que
    Que_nChrom_iLen_dic = {}
    for sID in sID_list:
      Que_chrsize_fNAME = chrsize_iPATH+'chrsize_'+sID+'.txt' # hard coding for chrsize files - need to be updated
      fpin = open(Que_chrsize_fNAME,'r')  
      for line in fpin:
        part = line.strip().split("\t")
        Que_nChrom = sID+"_"+part[0]
        iLen = part[1]
        Que_nChrom_iLen_dic.setdefault(Que_nChrom,iLen)
      fpin.close()

    # Writing_chrsize_chrname_chrcolor
    fpout_chrsize = open(BuscoCircosInput_byChr_oPATH+nLinkType+"ChrSize_"+ChrID+".tsv",'w')
    fpout_chrname = open(BuscoCircosInput_byChr_oPATH+nLinkType+"ChrName_"+ChrID+".tsv",'w')
    fpout_chrcolor= open(BuscoCircosInput_byChr_oPATH+nLinkType+"ChrColor_"+ChrID+".tsv",'w')

    nHeader_chrsize = "chrom\tchromStart\tchromEnd\tname\tgieStain\n"
    fpout_chrsize.write(nHeader_chrsize)

    nHeader_chrname = "chr\tlocus\tvalue\n"
    fpout_chrname.write(nHeader_chrname)

    nHeader_chrcolor = "sID\tcolor\n"
    fpout_chrcolor.write(nHeader_chrcolor)

    iCNT_chrline = 0
    for Ref_chr in Ref_Chr_list:
      iCNT_chrline += 1
      chriLen = Ref_nChrom_iLen_dic[Ref_chr]
      tmpline = Ref_chr+'\t'+'1'+'\t'+chriLen+'\t'+ Ref_id+'\t'+Ref_id+"\n"
      fpout_chrsize.write(tmpline)

      pos_in_circosplot = str(math.trunc(float(chriLen)/2))
      #pos_in_circosplot = str(chriLen)
      tmpline = Ref_chr +'\t'+ pos_in_circosplot  +'\t'+ Ref_chr + '\n'
      fpout_chrname.write(tmpline)

      chrID = str(iCNT_chrline)
      tmpline = chrID + '\t' + "1" + '\n'
      fpout_chrcolor.write(tmpline)

    #Que_Chr_list.reverse()
    for Que_chr in Que_Chr_list:
      sID = Que_chr.split("_")[0]
      nColor_sID = str(sID_list.index(sID)+2)
      try:
        iCNT_chrline += 1
        chriLen = Que_nChrom_iLen_dic[Que_chr]
        tmpline = Que_chr+'\t'+'1'+'\t'+chriLen+'\t'+ sID+'\t'+sID+"\n"
        fpout_chrsize.write(tmpline)

        pos_in_circosplot = str(math.trunc(float(chriLen)/2))
        #pos_in_circosplot = str(chriLen)
        tmpline = Que_chr +'\t'+ pos_in_circosplot  +'\t'+ Que_chr + '\n'
        fpout_chrname.write(tmpline)
        
        chrID = str(iCNT_chrline)
        tmpline = chrID + '\t' + nColor_sID + '\n'
        fpout_chrcolor.write(tmpline)
      
      except:
        print("Error-FASTA_BUSCO_ID_unmatched:",sID, Que_chr, Que_nChrom_iLen_dic.keys())
        sys.exit()
        fpout_chrsize.close()
        fpout_chrname.close()
        fpout_chrcolor.close()
    fpout_chrsize.close()
    fpout_chrname.close()
    fpout_chrcolor.close()
Generating_chrsize_chrname_chrcolor_linkcolor_from_link("Full")
Generating_chrsize_chrname_chrcolor_linkcolor_from_link("Core")


# Changing color codes ChrColor and LinkColor
def Changing_LinkColor_code(nTypeLink):# ex) nLinkType = Full
  flist = glob.glob(BuscoCircosInput_byChr_oPATH+nTypeLink+"LinkColor_*.tsv")
  for fNAME in flist:
    #nChrID = os.path.basename(fNAME).split("_")[1].split(".")[0]
    nID_iColorCode_dic = {}
    iColorCode_list = []
    nID_list = []
    fpin = open(fNAME,'r')
    nHeader = fpin.readline()
    for line in fpin:
      part = line.strip('\n').split('\t')
      nID= part[0]
      iColorCode = int(part[1])
      nID_list.append(nID)
      nID_iColorCode_dic.setdefault(nID,iColorCode)
      if not iColorCode == 0:
        if not iColorCode in iColorCode_list:
          iColorCode_list.append(iColorCode)
    fpin.close()

    EvenCode_list = iColorCode_list[math.ceil(len(iColorCode_list)/2.0):]
    iCNT_EvenCode = -1
    fpout = open(fNAME,'w')
    fpout.write(nHeader)
    for nID in nID_list:
      iColorCode = nID_iColorCode_dic[nID]
      if not iColorCode == math.ceil(iColorCode/2.0)*2: # Odd code
        nColorCode = str(int(math.ceil(iColorCode/2.0)))
      else:
        nColorCode = str(EvenCode_list[int((iColorCode/2)-1)])
      tmpline = nID+'\t'+nColorCode+'\n'
      fpout.write(tmpline)
    fpout.close()
Changing_LinkColor_code("Full")
Changing_LinkColor_code("Core")


def Changing_ChrColor_code(nLinkType):
  flist = glob.glob(BuscoCircosInput_byChr_oPATH+nLinkType+"ChrColor_*.tsv")
  for fNAME in flist:
    #nChrID = os.path.basename(fNAME).split("_")[1].split(".")[0]
    nID_iColorCode_dic = {}
    iColorCode_list = []
    nID_list = []
    fpin = open(fNAME,'r')
    nHeader = fpin.readline()
    for line in fpin:
      part = line.strip('\n').split('\t')
      nID= part[0]
      iColorCode = int(part[1])
      nID_list.append(nID)
      nID_iColorCode_dic.setdefault(nID,iColorCode)
      if not iColorCode == 0:
        if not iColorCode in iColorCode_list:
          iColorCode_list.append(iColorCode)
    fpin.close()
    EvenCode_list = iColorCode_list[1:math.trunc(len(iColorCode_list)/2.0)+1]
    OddCode_list = iColorCode_list[math.trunc(len(iColorCode_list)/2.0)+1:]
    iCNT_EvenCode = -1
    fpout = open(fNAME,'w')
    fpout.write(nHeader)
    for nID in nID_list:
      iColorCode = nID_iColorCode_dic[nID]
      if iColorCode == 1:
        nColorCode = str(1)
      else:
        if iColorCode == math.ceil(iColorCode/2.0)*2: # Even code
          nColorCode = str(EvenCode_list[int(iColorCode/2.0)-1])
        else: # Odd code >= 3
          nColorCode = str(OddCode_list[int(iColorCode/2.0)-1])
      tmpline = nID+'\t'+nColorCode+'\n'
      fpout.write(tmpline)
    fpout.close()
Changing_ChrColor_code("Full")
Changing_ChrColor_code("Core")






'''
# Generating_chrsize_chrname_chrcolor_linkcolor_from_link
link_flist = glob.glob(BuscoCircosInput_oPATH+"Link_*.tsv")
for fNAME_link in link_flist:
  # reading link files and setting list of chr of each sID
  # Writing_link color by chr
  sID = os.path.basename(fNAME_link).split("_")[1].split(".")[0]

  LinkColor_fNAME = BuscoCircosInput_oPATH+"LinkColor_"+sID+'.tsv'
  fpout_LinkColor = open(LinkColor_fNAME,'w')
  LinkColor_nHeader = "link_ID\tcolor\n"
  fpout_LinkColor.write(LinkColor_nHeader)
  
  fpin = open(fNAME_link,'r')
  fpin.readline()
  Ref_Chr_list = []
  Que_Chr_list = []
  iCNT_linkID = 0
  iColorCode_Link = 0
  for line in fpin:
    iCNT_linkID += 1
    nCNT_linkID = str(iCNT_linkID)
    part = line.strip('\n').split("\t")
    Ref_chr = part[0]
    if not Ref_chr in Ref_Chr_list:
      iColorCode_Link += 1
      nColorCode_Link = str(iColorCode_Link)
      Ref_Chr_list.append(Ref_chr)
    Que_chr = part[3]
    if not Que_chr in Que_Chr_list:
      Que_Chr_list.append(Que_chr)
    tmpline = nCNT_linkID + '\t' + nColorCode_Link + '\n'
    fpout_LinkColor.write(tmpline)
  fpin.close()
  fpout_LinkColor.close()
  
  # Reading_chrsize_Ref
  Ref_id = nRefAsmID.split("_")[1]
  Ref_chrsize_fNAME = chrsize_iPATH+'chrsize_'+Ref_id+'.txt' # hard coding for chrsize files - need to be updated
  fpin = open(Ref_chrsize_fNAME,'r') 
  Ref_nChrom_iLen_dic = {}
  for line in fpin:
    part = line.strip().split("\t")
    Ref_nChrom_iLen_dic.setdefault(part[0],part[1])
  fpin.close()

  # Reading_chrsize_Que
  Que_chrsize_fNAME = chrsize_iPATH+'chrsize_'+sID+'.txt' # hard coding for chrsize files - need to be updated
  fpin = open(Que_chrsize_fNAME,'r')  
  Que_nChrom_iLen_dic = {}
  for line in fpin:
    part = line.strip().split("\t")
    Que_nChrom_iLen_dic.setdefault(part[0],part[1])
  fpin.close()

  # Writing_chrsize_chrname_chrcolor
  fpout_chrsize = open(BuscoCircosInput_oPATH+"ChrSize_"+sID+".tsv",'w')
  fpout_chrname = open(BuscoCircosInput_oPATH+"ChrName_"+sID+".tsv",'w')
  fpout_chrcolor= open(BuscoCircosInput_oPATH+"ChrColor_"+sID+".tsv",'w')

  nHeader_chrsize = "chrom\tchromStart\tchromEnd\tname\tgieStain\n"
  fpout_chrsize.write(nHeader_chrsize)

  nHeader_chrname = "chr\tlocus\tvalue\n"
  fpout_chrname.write(nHeader_chrname)

  nHeader_chrcolor = "sID\tcolor\n"
  fpout_chrcolor.write(nHeader_chrcolor)

  iCNT_chrline = 0
  for Ref_chr in Ref_Chr_list:
    iCNT_chrline += 1
    chriLen = Ref_nChrom_iLen_dic[Ref_chr]
    tmpline = Ref_chr+'\t'+'1'+'\t'+chriLen+'\t'+ Ref_id+'\t'+Ref_id+"\n"
    fpout_chrsize.write(tmpline)

    pos_in_circosplot = str(math.trunc(float(chriLen)/2))
    tmpline = Ref_chr +'\t'+ pos_in_circosplot  +'\t'+ Ref_chr + '\n'
    fpout_chrname.write(tmpline)

    chrID = str(iCNT_chrline)
    tmpline = chrID + '\t' + "1" + '\n'
    fpout_chrcolor.write(tmpline)

  Que_Chr_list.reverse()
  for Que_chr in Que_Chr_list:
    try:
      iCNT_chrline += 1
      chriLen = Que_nChrom_iLen_dic[Que_chr]
      tmpline = Que_chr+'\t'+'1'+'\t'+chriLen+'\t'+ sID+'\t'+sID+"\n"
      fpout_chrsize.write(tmpline)

      pos_in_circosplot = str(math.trunc(float(chriLen)/2))
      tmpline = Que_chr +'\t'+ pos_in_circosplot  +'\t'+ Que_chr + '\n'
      fpout_chrname.write(tmpline)
      
      chrID = str(iCNT_chrline)
      tmpline = chrID + '\t' + "2" + '\n'
      fpout_chrcolor.write(tmpline)
    
    except:
      print("Error-FASTA_BUSCO_ID_unmatched:",sID, Que_chr)
      sys.exit()
      fpout_chrsize.close()
      fpout_chrname.close()
      fpout_chrcolor.close()
  fpout_chrsize.close()
  fpout_chrname.close()
  fpout_chrcolor.close()
'''

  
