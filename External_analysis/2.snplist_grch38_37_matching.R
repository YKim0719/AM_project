
aa=fread("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/Insomnia/sumstat_pheno_all.txt",header=T,stringsAsFactors=FALSE)
aa=as.data.frame(aa)
chr=1
dat2=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/prs_external/EA_ldpred/info2/chr_",chr,"_ver2.txt",sep=""),header=T,stringsAsFactors=FALSE)
dat2=as.data.frame(dat2)
if(sum(colnames(dat2)=="GRCH38_pos")>=2){
  dat2=dat2[,-ncol(dat2)]
}
bim=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/new_LD/cognitive_all/chr.",chr,".bim",sep=""),header=F,stringsAsFactors=FALSE)
bim=as.data.frame(bim)
colnames(bim)=c("chr","SNPID_GRCH38","MISS","GRCH38_pos","REF","ALT")
bim=bim[,c(2,4)]
dat_merged_ref=merge(dat2,bim,by="GRCH38_pos")

for(chr in 2:22){
  dat2=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/prs_external/EA_ldpred/info2/chr_",chr,"_ver2.txt",sep=""),header=T,stringsAsFactors=FALSE)
  dat2=as.data.frame(dat2)
  if(sum(colnames(dat2)=="GRCH38_pos")>=2){
    dat2=dat2[,-ncol(dat2)]
  }
  bim=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/new_LD/cognitive_all/chr.",chr,".bim",sep=""),header=F,stringsAsFactors=FALSE)
  bim=as.data.frame(bim)
  colnames(bim)=c("chr","SNPID_GRCH38","MISS","GRCH38_pos","REF","ALT")
  bim=bim[,c(2,4)]
  dat_merged_ref2=merge(dat2,bim,by="GRCH38_pos")
  dat_merged_ref=rbind(dat_merged_ref,dat_merged_ref2)
}

head(dat_merged_ref)
dat_merged_ref2=dat_merged_ref[,c(2,3,1,13,16)]
fwrite(dat_merged_ref2,"/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/matched_snp_info.txt",row.names=F,col.names=T,quote=F,sep="\t")

dat_matched=dat[match(dat_merged_ref2$rsid,dat$RSID),]


cd /pl/active/KellerLab/Yongkang/UKBhap/prs_external/summary_statistics/gscan
gunzip SmokingInitiation.WithoutUKB.txt.gz
gunzip SmokingCessation.WithoutUKB.txt.gz
gunzip DrinksPerWeek.WithoutUKB.txt.gz
gunzip CigarettesPerDay.WithoutUKB.txt.gz
gunzip AgeOfInitiation.WithoutUKB.txt.gz
