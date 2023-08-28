#ml anaconda
#conda activate python3conda
#conda create -n python3conda python==3.6.8
#conda install joblib 
#conda install matplotlib.path
#conda install paratext pandas
#python -m pip install pandas
#chr=str(1);idx=str(1)
#phenonames=c("AgeOfInitiation_gscan","CigarettesPerDay_gscan","DrinksPerWeek_gscan","SmokingCessation_gscan","SmokingInitiation_gscan")
import numpy as np
import pandas as pd
import pickle
import os
import sys
#phenoname="AgeOfInitiation_gscan"
chr = sys.argv[1]
idx = sys.argv[2]
dirname = sys.argv[3]
#chr, idx , phenoname = input().split()

#chr=chr.apply("str")
#idx=idx.apply("str")
#phenoname=phenoname.apply("str")
#chr=chr.astype(str)
#idx=idx.astype(str)
#phenoname=phenoname.astype(str)
chr2=int(chr)

ids1 = np.loadtxt("/pl/active/IBG/UKB_TOPMed_imputed_pending/vcf_GT/hapmap3_haps/split_chr_vcf/UKB_chr"+chr+".idx"+idx+".indep.haps",skiprows=0,dtype="str",usecols=(0,1,2,3,4,5),delimiter="\t")
id = np.loadtxt("/pl/active/IBG/UKB_TOPMed_imputed_pending/vcf_GT/hapmap3_haps/split_chr_vcf/UKB_chr"+chr+".idx"+idx+".indep.haps",skiprows=0,dtype="bool_",delimiter="\t",converters = {0: str, 1: str, 2: str, 3: str, 4: str, 5: str, 6: str, 7: str, 8: str})

  
#ids1 = np.loadtxt("/pl/active/KellerLab/Yongkang/UKBhap/Meng/haps/indep_samp_hap/UKB_chr"+chr+".idx"+idx+".indep.haps",skiprows=1,dtype="str",usecols=(0,1,2,3,4,5),delimiter=",")
#id = np.loadtxt("/pl/active/KellerLab/Yongkang/UKBhap/Meng/haps/indep_samp_hap/UKB_chr"+chr+".idx"+idx+".indep.haps",skiprows=1,dtype="bool_",delimiter=",",converters = {0: str, 1: str
#, 2: str, 3: str, 4: str, 5: str, 6: str, 7: str, 8: str})
pgs_score_inf=np.loadtxt(dirname+"/ldpred_results.csv",delimiter=",",dtype="float",skiprows=1,usecols=11)
#pgs_score_auto=np.loadtxt(dirname+"/ldpred_results.csv",delimiter=",",dtype="float",skiprows=1,usecols=12)
ids2=np.loadtxt(dirname+"/ldpred_results.csv",delimiter=",",dtype="<U23",skiprows=1,usecols=0)
ref2=np.loadtxt(dirname+"/ldpred_results.csv",delimiter=",",dtype="str",skiprows=1,usecols=4)

iid = np.loadtxt("/pl/active/IBG/UKB_TOPMed_imputed_pending/vcf_GT/hapmap3_haps/split_chr_vcf/IID_list.txt",dtype="float",comments='!')
from typing import Hashable, List
def match_list(a: List[Hashable], b: List[Hashable]) -> List[int]:
  return [b.index(x) if x in b else None for x in a]

from typing import Hashable, List
def include(a: List[Hashable], b: List[Hashable]) -> List[bool]:
  return [True if x in b else False for x in a]

idx1=np.ndarray.tolist(ids1[:,2])
idx2=np.ndarray.tolist(ids2)
aa=include(idx1,idx2)
ids1=ids1[aa,:]
hap_dat=id[aa,:]
index=np.r_[0:9]
hap_dat=np.delete(hap_dat, index,1)
ref1=ids1[:,3]
idx1=np.ndarray.tolist(ids1[:,2])

bb=match_list(idx1,idx2)
bb=np.array(bb, dtype=np.float)
bb=bb.astype('int')
pgs_score_inf2=pgs_score_inf[bb]
#pgs_score_auto2=pgs_score_auto[bb]
ids2=ids2[bb]
idx2=np.ndarray.tolist(ids2)
ref2=ref2[bb]
if sum(ref1!=ref2)!=0 :
  pgs_score_inf2[ref1!=ref2]=-pgs_score_inf2[ref1!=ref2]
#  pgs_score_auto2[ref1!=ref2]=-pgs_score_auto2[ref1!=ref2]

result_inf=np.matmul(pgs_score_inf2,hap_dat)
#åresult_auto=np.matmul(pgs_score_auto2,hap_dat)
#result=np.dot(hap_dat[:,0],pgs_score2)
odd=np.linspace(start=0, stop=len(result_inf)-2, num=len(result_inf)/2)
even=np.linspace(start=1, stop=len(result_inf)-1, num=len(result_inf)/2)
odd=odd.astype("int")
even=even.astype("int")
iid=iid.astype("int")
results=pd.DataFrame(
  {"IID": iid,
  "pgs_inf1": result_inf[odd],
  "pgs_inf2": result_inf[even],
#  "pgs_auto1": result_auto[odd],
#  "pgs_auto2": result_auto[even]
  }
)
#"
  
dirName=dirname+"/individual_hpgs"
if not os.path.exists(dirName):
    os.makedirs(dirName)

#os.makedirs("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/"+phenoname+"/whole_variant/individual_hpgs2")
np.savetxt(dirName+"/chr"+chr+"."+idx+"_ver2.csv",results,delimiter=",")

