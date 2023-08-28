args <- commandArgs(trailingOnly =TRUE)
chr=as.numeric(args[1])
idx=as.numeric(args[2])
library(bigsnpr)

library(data.table)

if(chr<=12){
  dat=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/haps/indep_samp_hap/UKB_chr",chr,".idx",idx,".indep.haps",sep=""),header=T,stringsAsFactors=FALSE)  
}else{
  dat=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/haps/indep_samp_hap/UKB_chr",chr,".idx",idx,".indep.haps",sep=""),header=F,stringsAsFactors=FALSE)  
  
}
dir.create("/rc_scratch/yoki5348/indep_samp_hap/")
dat=as.data.frame(dat)

pgs=fread("/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/height/ldpred_results.csv",header=T,stringsAsFactors=FALSE)
pgs=as.data.frame(pgs)

dat=dat[which(dat$V3%in%pgs$ID),]

fwrite(dat,paste("/rc_scratch/yoki5348/indep_samp_hap/UKB_chr",chr,".idx",idx,".indep.hapmap3.haps",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
