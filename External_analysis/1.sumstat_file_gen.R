# Unify format of summary statistics from the biobanks or references to the format which can be used in LDSC
library(data.table)
dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/intelligence/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt",header=T,stringsAsFactors=FALSE)
dat=as.data.frame(dat)
dat$A1=toupper(dat$A1) # To change whole allele information to capital cases.
dat$A2=toupper(dat$A2)
colnames(dat)[c(1,3,4)]=c("rsid","chr","pos")

snpinfo=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/matched_snp_info.txt",header=T,stringsAsFactors=FALSE)
# Since Haps files were imputed by GRCH38 reference panels, we changed GRCH37 summary information into GRCH37 format.
snpinfo=as.data.frame(snpinfo)
dat=merge(dat,snpinfo,by=c("rsid","chr"))
dat2=fread("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/smoke/sumstat_pheno_all.txt",header=T,stringsAsFactors=FALSE)
dat2=as.data.frame(dat2)
dir.create("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/intelligence")
dat3=dat[,c(17,2,16,5,6,9,10,11,13,12)]
colnames(dat3)=c(colnames(dat2),"N") 
fwrite(dat3,"/pl/active/KellerLab/Yongkang/PGS_web/ldsc/intelligence/sumstat_pheno_all.txt",row.names=F,col.names=T,quote=F,sep="\t")

dat=fread("/pl/active/KellerLab/Yongkang/UKBhap/prs_external/summary_statistics/GWAS_CP_all.txt",header=T,stringsAsFactors=FALSE)
dat=as.data.frame(dat)
colnames(dat)[c(1,2,3)]=c("rsid","chr","pos")
snpinfo=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/matched_snp_info.txt",header=T,stringsAsFactors=FALSE)
snpinfo=as.data.frame(snpinfo)
dat=merge(dat,snpinfo,by=c("rsid","chr"))
dat2=fread("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/smoke/sumstat_pheno_all.txt",header=T,stringsAsFactors=FALSE)
dat2=as.data.frame(dat2)
dat3=dat[,c(12,2,11,4,5,7,8,9)]
dat3$info=0.95 #Since this summary statistic does not contain INFO score, put 0.95 for whole SNPs.
colnames(dat3)=colnames(dat2)
fwrite(dat3,"/pl/active/KellerLab/Yongkang/PGS_web/ldsc/cognitive/sumstat_pheno_all.txt",row.names=F,col.names=T,quote=F,sep="\t")
