fpin = open("reordered_Core_Matrix_BUSCO_sorted.csv",'r')
header= fpin.readline().strip().split(",")
sID_list = header[1:]

ChrID_sID_otherChrID_dic = {}
for line in fpin:
  part = line.strip().split(",")
  ChrID = part[1].split("&")[0]
  ChrID_sID_otherChrID_dic.setdefault(ChrID,{})
  iCNT_sID = len(sID_list)+1
  for i in range(2,iCNT_sID):
    sID = sID_list[i-1]
    ChrID_sID_otherChrID_dic[ChrID].setdefault(sID, {})
    otherChrID = part[i].split("&")[0]
    ChrID_sID_otherChrID_dic[ChrID][sID].setdefault(otherChrID, '')
fpin.close()

fpout = open("summary_reordered_Core_Matrix_BUSCO_sorted.csv",'w')
header = '\t'.join(sID_list)+'\t'+'\t'.join(sID_list[1:])+'\n'
fpout.write(header)
for ChrID in ChrID_sID_otherChrID_dic.keys():
  tmpline = ChrID
  for sID in sID_list[1:]:
    tmpline += '\t' + str(len(ChrID_sID_otherChrID_dic[ChrID][sID].keys()))
  for sID in sID_list[1:]:
    tmp_list = []
    for otherChrID in ChrID_sID_otherChrID_dic[ChrID][sID].keys():
      tmp_list.append(otherChrID)
    tmpline += '\t'+','.join(tmp_list)
  fpout.write(tmpline+'\n')
fpout.close()
    
  
