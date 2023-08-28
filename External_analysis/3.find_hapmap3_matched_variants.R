###Get Hapamp3 mapped SNPs information from NCBI######
for chr in {1..22}
do
wget "ftp.ncbi.nlm.nih.gov/hapmap/frequencies/2010-05_phaseIII/allele_freqs_chr"$chr"_CEU_phase3.3_nr.b36_fwd.txt.gz"
done

for chr in {1..22}
do
gunzip "allele_freqs_chr"$chr"_CEU_phase3.3_nr.b36_fwd.txt.gz"
done
more "allele_freqs_chr"$chr"_CEU_phase3.3_nr.b36_fwd.txt"

####################Since Meng's imputed results are built from GRCH38 build, we had to match them to b36 position which is used for hapmap3 NCBI DB#####
chr=1
library(data.table)
dat=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/matched_snp_info.txt",header=T,stringsAsFactors = F)
dat=as.data.frame(dat2)

for(chr in 1:22){
  dat=read.table(paste("/pl/active/KellerLab/Yongkang/hapmap3/allele_freqs_chr",chr,"_CEU_phase3.3_nr.b36_fwd.txt",sep=""),header=F,skip=1,stringsAsFactors=FALSE,sep=" ")
  dat=dat[which(dat$V1%in%dat2$rsid),]  
  dat=cbind(dat[,c(1,3,11,14,12)],dat2[match(dat$V1,dat2$rsid),])
  colnames(dat)[1:5]=c("RSID","B35_POS","A1","A2","EAF_A1")
  write.table(dat,paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/hapmap3_chr",chr,".txt",sep=""),row.names=F,col.names=T,quote=F)
}

chr=1
dat2=fread(paste("/pl/active/KellerLab/Yongkang/hapmap3/allele_freqs_chr",chr,"_CEU_phase3.3_nr.b36_fwd.txt",sep=""),header=F,skip=1,stringsAsFactors=FALSE,sep=" ")
dat2=as.data.frame(dat2)
aa=paste(dat2$V2,":",dat2$V3,"-",dat2$V3,sep="")
for(chr in 2:22){
  dat3=fread(paste("/pl/active/KellerLab/Yongkang/hapmap3/allele_freqs_chr",chr,"_CEU_phase3.3_nr.b36_fwd.txt",sep=""),header=F,skip=1,stringsAsFactors=FALSE,sep=" ")
  dat3=as.data.frame(dat3)
  aa=c(aa,paste(dat3$V2,":",dat3$V3,"-",dat3$V3,sep=""))
  dat2=rbind(dat2,dat3)
}
hapmap_all=dat2
write.table(aa,"/projects/yoki5348/all_hapmap3.txt",row.names=F,col.names=F,quote=F)
write.table(hapmap_all,"/projects/yoki5348/hapmap3_info_all.txt",row.names=F,col.names=T,quote=F)

aa=fread("/projects/yoki5348/all_hapmap3.txt",header=F,stringsAsFactors=FALSE)
aa=as.data.frame(aa)[,1]
hapmap_all=fread("/projects/yoki5348/hapmap3_info_all.txt",header=T,stringsAsFactors=FALSE)
hapmap_all=as.data.frame(hapmap_all)
dat2=fread("/projects/yoki5348/hglft_genome_51b3_65c9c0.bed",header=F,stringsAsFactors=FALSE)
dat2=as.data.frame(dat2)
dat3=fread("/projects/yoki5348/hglft_genome_51b3_65c9c0.err.txt",header=F,stringsAsFactors=FALSE)
dat3=dat3[seq(2,nrow(dat3),2),]
dat3=as.data.frame(dat3)
dat3=dat3[,1]
aa2=aa[-which(aa%in%dat3)]
b36_GRCH38=data.frame(b36=aa2,GRCH38=dat2[,1])
hapmap_all2=hapmap_all[-which(aa%in%dat3),]
b36_GRCH38=cbind(b36_GRCH38,hapmap_all2[,c(1,2,3,4,11,14)])
colnames(b36_GRCH38)[3:8]=c("rsid","chr","b36pos","strand","REF","ALT")
library(stringr)
pos=str_split(b36_GRCH38$GRCH38,"-")
pos2=c()
library(foreach)
pos2 <- foreach(idx=1:length(pos), .combine='c') %do% as.numeric(pos[[idx]][2])
b36_GRCH38=cbind(b36_GRCH38,pos2)
colnames(b36_GRCH38)[9]="GRCH38pos"
write.table(b36_GRCH38,"/pl/active/KellerLab/Yongkang/hapmap3/GRCH38_matched.txt",row.names=F,col.names=T,quote=F,sep="\t")
b36_GRCH38=fread("/pl/active/KellerLab/Yongkang/hapmap3/GRCH38_matched.txt",header=T,stringsAsFactors=FALSE)
b36_GRCH38=as.data.frame(b36_GRCH38)

chr=1
hapmap_all3=b36_GRCH38[which(b36_GRCH38$chr==paste("chr",chr,sep="")),]
bb=paste(hapmap_all3$chr,":",hapmap_all3$GRCH38pos,":",hapmap_all3$REF,":",hapmap_all3$ALT,sep="")
bb2=paste(hapmap_all3$chr,":",hapmap_all3$GRCH38pos,":",hapmap_all3$ALT,":",hapmap_all3$REF,sep="")
###############################################################
neg_strand=rep(NA,length(hapmap_all3$REF))###Strand match
neg_strand[hapmap_all3$REF=="A"]="T"
neg_strand[hapmap_all3$REF=="T"]="A"
neg_strand[hapmap_all3$REF=="G"]="C"
neg_strand[hapmap_all3$REF=="C"]="G"

neg_strand2=rep(NA,length(hapmap_all3$ALT))
neg_strand2[hapmap_all3$ALT=="A"]="T"
neg_strand2[hapmap_all3$ALT=="T"]="A"
neg_strand2[hapmap_all3$ALT=="G"]="C"
neg_strand2[hapmap_all3$ALT=="C"]="G"
bb3=paste(hapmap_all3$chr,":",hapmap_all3$GRCH38pos,":",neg_strand,":",neg_strand2,sep="")
bb4=paste(hapmap_all3$chr,":",hapmap_all3$GRCH38pos,":",neg_strand2,":",neg_strand,sep="")

pgen=fread(paste("/pl/active/IBG/UKB_TOPMed_imputed_pending/vcf_GT/chr",chr,"_v3.pvar",sep=""),header=T,stringsAsFactors=FALSE)
pgen=as.data.frame(pgen)
pgen=pgen[which(pgen$REF%in%c("A","T","G","C")),]
pgen=pgen[which(pgen$ALT%in%c("A","T","G","C")),]#Remove structural variations

results=rep(NA,length(bb))
results[which(bb3%in%pgen$ID)]=bb3[which(bb3%in%pgen$ID)]
results[which(bb4%in%pgen$ID)]=bb4[which(bb4%in%pgen$ID)]
results[which(bb%in%pgen$ID)]=bb[which(bb%in%pgen$ID)]
results[which(bb2%in%pgen$ID)]=bb2[which(bb2%in%pgen$ID)]
hapmap_all3$MENG_ID=results


write.table(results,paste("/pl/active/KellerLab/Yongkang/hapmap3/chr",chr,"GRCH38_matched_MENG.txt",sep=""),row.names=F,col.names=F,quote=F)

for(chr in 2:22){
  hapmap_all4=b36_GRCH38[which(b36_GRCH38$chr==paste("chr",chr,sep="")),]
  bb=paste(hapmap_all4$chr,":",hapmap_all4$GRCH38pos,":",hapmap_all4$REF,":",hapmap_all4$ALT,sep="")
  bb2=paste(hapmap_all4$chr,":",hapmap_all4$GRCH38pos,":",hapmap_all4$ALT,":",hapmap_all4$REF,sep="")
  
  neg_strand=rep(NA,length(hapmap_all4$REF))
  neg_strand[hapmap_all4$REF=="A"]="T"
  neg_strand[hapmap_all4$REF=="T"]="A"
  neg_strand[hapmap_all4$REF=="G"]="C"
  neg_strand[hapmap_all4$REF=="C"]="G"
  
  neg_strand2=rep(NA,length(hapmap_all4$ALT))
  neg_strand2[hapmap_all4$ALT=="A"]="T"
  neg_strand2[hapmap_all4$ALT=="T"]="A"
  neg_strand2[hapmap_all4$ALT=="G"]="C"
  neg_strand2[hapmap_all4$ALT=="C"]="G"
  bb3=paste(hapmap_all4$chr,":",hapmap_all4$GRCH38pos,":",neg_strand,":",neg_strand2,sep="")
  bb4=paste(hapmap_all4$chr,":",hapmap_all4$GRCH38pos,":",neg_strand2,":",neg_strand,sep="")
  
  pgen=fread(paste("/pl/active/IBG/UKB_TOPMed_imputed_pending/vcf_GT/chr",chr,"_v3.pvar",sep=""),header=T,stringsAsFactors=FALSE)
  pgen=as.data.frame(pgen)
  pgen=pgen[which(pgen$REF%in%c("A","T","G","C")),]
  pgen=pgen[which(pgen$ALT%in%c("A","T","G","C")),]
  
  results=rep(NA,length(bb))
  results[which(bb3%in%pgen$ID)]=bb3[which(bb3%in%pgen$ID)]
  results[which(bb4%in%pgen$ID)]=bb4[which(bb4%in%pgen$ID)]
  results[which(bb%in%pgen$ID)]=bb[which(bb%in%pgen$ID)]
  results[which(bb2%in%pgen$ID)]=bb2[which(bb2%in%pgen$ID)]
  
  write.table(results,paste("/pl/active/KellerLab/Yongkang/hapmap3/chr",chr,"GRCH38_matched_MENG.txt",sep=""),row.names=F,col.names=F,quote=F)
  hapmap_all4$MENG_ID=results
  hapmap_all3=rbind(hapmap_all3,hapmap_all4)
  
}

aa=c()
for(chr in 1:22){
  results=fread(paste("/pl/active/KellerLab/Yongkang/hapmap3/chr",chr,"GRCH38_matched_MENG.txt",sep=""),header=F,stringsAsFactors=FALSE)
  results=as.data.frame(results)
  aa=c(aa,as.character(results[,1]))
}
hapmap_all3$MENG_ID=as.character(hapmap_all3$MENG_ID)
hapmap_all3$MENG_ID=aa

fwrite(hapmap_all3,"/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",row.names=F,col.names=T,quote=F,sep="\t")

sum(pgen$ID%in%bb)+sum(pgen$ID%in%bb2)+
sum(pgen$ID%in%bb3)+
sum(pgen$ID%in%bb4)
snplist=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/chr_22_snplist.txt",header=F,stringsAsFactors=FALSE)
pgen=fread("/pl/active/IBG/UKB_TOPMed_imputed_pending/vcf_GT/chr22_v3.pvar",header=T,stringsAsFactors=FALSE)
pgen=as.data.frame(pgen)
pgen=pgen[which(pgen$REF%in%c("A","T","G","C")),]
pgen=pgen[which(pgen$ALT%in%c("A","T","G","C")),]
pgen_pos=paste("chr",pgen[,1],":",pgen[,2],"-",pgen[,2],sep="")
sum(pgen_pos%in%bb)
sum(pgen_pos%in%bb2)
sum(pgen_pos%in%bb3)
sum(pgen_pos%in%bb4)
pgen2=pgen[which(pgen_pos%in%bb),]


pgen2[pgen2[,2]==16802212,]
library(data.table)
chr=1
dat2=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/hapmap3_chr",chr,".txt",sep=""),header=T,stringsAsFactors=FALSE)
dat=as.data.frame(dat)
dim(dat)
a=nrow(dat)
for(chr in 2:22){
  dat=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/hapmap3_chr",chr,".txt",sep=""),header=T,stringsAsFactors=FALSE)
  dat=as.data.frame(dat)
  a=c(a,nrow(dat))
  
}
sum(dat$V1%in%dat2$rsid)
dat3=dat2[dat2$chr==1,]

for(chr in 1:22){
  system(paste("~/plink2 --pfile /pl/active/IBG/UKB_TOPMed_imputed_pending/vcf_GT/chr",chr,"_v3 --extract /pl/active/KellerLab/Yongkang/hapmap3/chr",chr,"GRCH38_matched_MENG.txt ",
               "--keep /pl/active/KellerLab/Yongkang/UKBhap/prs_external/LD_files/individual.txt ",
               "--make-bed --out /pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,sep=""))
}



###########################Perform LDPRED2 ###############################
dir.create("/pl/active/KellerLab/Yongkang/hapmap3/ldfile")
NCORES <- 10
dir.create("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data")
tmp <- tempfile(tmpdir = "/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)
for(chr in 1:22){
  try(snp_readBed(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".bed",sep="")))
  obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  # perform SNP matching
  genotype <- obj.bigSNP$genotypes
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  # calculate LD
  
  # Extract SNPs that are included in the chromosome
  ind.chr2 <- 1:nrow(map)
  # Calculate the LD
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = NCORES,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  ld <- Matrix::colSums(corr0^2)
  map=cbind(map,ld)
  dir.create("/pl/active/KellerLab/Yongkang/hapmap3/ldfile")
  dir.create("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data")
  save(list="corr0",file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
  
  fwrite(map,paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/new_LD/indep_UKB/ldfile/tmp-data/mapfile_chr",chr,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  
}



######################################perform LDpred - Bipolar disorder ###############################


dat=fread("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/bipolar/bipolar1/pgc-bip2021-all.vcf.tsv",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("ID","A1","A2","BETA","SE","NEFFDIV2")]
colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

chr=1
mapfile=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/new_LD/indep_UKB/ldfile/tmp-data/mapfile_chr",chr,".txt",sep=""),header=T,stringsAsFactors=FALSE)
mapfile=as.data.frame(mapfile)
mapfile$"_NUM_ID_"=1:nrow(mapfile)
mapfile$"NUM_each_chr"=1:nrow(mapfile)
for(chr in 2:22){
  mapfile2=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/new_LD/indep_UKB/ldfile/tmp-data/mapfile_chr",chr,".txt",sep=""),header=T,stringsAsFactors=FALSE)
  mapfile2=as.data.frame(mapfile2)
  mapfile2$"_NUM_ID_"=nrow(mapfile)+(1:nrow(mapfile2))
  mapfile2$"NUM_each_chr"=1:nrow(mapfile2)
  mapfile=rbind(mapfile,mapfile2)
}
colnames(mapfile)=c("CHR","ID","POS","A0","A1","ld","_NUM_ID","NUM_each_chr")
fwrite(mapfile,"/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",row.names=F,col.names=T,quote=F,sep="\t")

mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]


NCORES <- 20
# Open a temporary file
tmp <- tempfile(tmpdir = "/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

chr=1
load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
dat_merged2=dat_merged[dat_merged$CHR==chr,]
corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
corr <- as_SFBM(corr0, tmp)
for(chr in 2:22){
  load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
  dat_merged2=dat_merged[dat_merged$CHR==chr,]
  corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
  corr$add_columns(corr0, nrow(corr))
}
save(list="corr",file="/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr_file_final_bipolar_file1.RData")  

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
chr=1
obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
G   <- obj.bigSNP$genotypes

beta_auto2=beta_auto[dat_merged$CHR==chr,]
df_beta2=df_beta[dat_merged$CHR==chr,]
pred_auto <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                         ind.col = df_beta2[["NUM_each_chr"]])

for(chr in 2:22){
  obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
  G   <- obj.bigSNP$genotypes
  
  beta_auto2=beta_auto[dat_merged$CHR==chr,]
  df_beta2=df_beta[dat_merged$CHR==chr,]
  pred_auto2 <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                           ind.col = df_beta2[["NUM_each_chr"]])
  pred_auto=pred_auto+pred_auto2
}

sc <- apply(pred_auto, 2, sd)
# Make TRUE FALSE vector for which converged
convergedAuto <- ((abs(sc - median(sc,na.rm=T)) < 3 * mad(sc,na.rm=T)) & !is.na(sc))
# Average the betas for the models that converged
final_beta_auto <- rowMeans(beta_auto[, convergedAuto])
dat_merged$beta_ldpred_auto=final_beta_auto
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/bipolar/bipolar1/ldpred_results.csv",row.names=F,quote=F)

############Bipolar2
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=fread("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/bipolar/bipolar2/daner_PGC_BIP32b_mds7a_0416a",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE","Neff")]
dat$OR=log(dat$OR)
colnames(dat)[4]="BETA"
colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile=as.data.frame(mapfile)

mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]


NCORES <- 10
# Open a temporary file
dir.create("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/bipolar_file2")
tmp <- tempfile(tmpdir = "/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/bipolar_file2")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

chr=1
load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
dat_merged2=dat_merged[dat_merged$CHR==chr,]
corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
corr <- as_SFBM(corr0, tmp)
for(chr in 2:22){
  load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
  dat_merged2=dat_merged[dat_merged$CHR==chr,]
  corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
  corr$add_columns(corr0, nrow(corr))
}
save(list="corr",file="/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/PTSD/corr_file_final.RData")  

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
chr=1
obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
G   <- obj.bigSNP$genotypes

beta_auto2=beta_auto[dat_merged$CHR==chr,]
df_beta2=df_beta[dat_merged$CHR==chr,]
pred_auto <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                         ind.col = df_beta2[["NUM_each_chr"]])

for(chr in 2:22){
  obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
  G   <- obj.bigSNP$genotypes
  
  beta_auto2=beta_auto[dat_merged$CHR==chr,]
  df_beta2=df_beta[dat_merged$CHR==chr,]
  pred_auto2 <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                            ind.col = df_beta2[["NUM_each_chr"]])
  pred_auto=pred_auto+pred_auto2
}

sc <- apply(pred_auto, 2, sd)
# Make TRUE FALSE vector for which converged
convergedAuto <- ((abs(sc - median(sc,na.rm=T)) < 3 * mad(sc,na.rm=T)) & !is.na(sc))
# Average the betas for the models that converged
final_beta_auto <- rowMeans(beta_auto[, convergedAuto])
dat_merged$beta_ldpred_auto=final_beta_auto

fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD1/ldpred_results.csv",row.names=F,quote=F)

######################Bipolar 3


library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=fread("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/bipolar/bipolar3/pgc.bip.full.2012-04.txt",header=F,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
colnames(dat)[1:7]=c("SNP","A1","POS","A1","A2","OR","SE")
dat=cbind(dat,400000)
colnames(dat)[ncol(dat)]="Neff"
dat=dat[,c("SNP","A1","A2","OR","SE","Neff")]
dat$OR=log(dat$OR)
colnames(dat)[4]="BETA"
colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile=as.data.frame(mapfile)

mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]


NCORES <- 10
# Open a temporary file
dir.create("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/bipolar_file3")
tmp <- tempfile(tmpdir = "/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/bipolar_file3")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

chr=1
load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
dat_merged2=dat_merged[dat_merged$CHR==chr,]
corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
corr <- as_SFBM(corr0, tmp)
for(chr in 2:22){
  load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
  dat_merged2=dat_merged[dat_merged$CHR==chr,]
  corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
  corr$add_columns(corr0, nrow(corr))
}
save(list="corr",file="/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/bipolar_file3/corr_file_final.RData")  

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
chr=1
obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
G   <- obj.bigSNP$genotypes

beta_auto2=beta_auto[dat_merged$CHR==chr,]
df_beta2=df_beta[dat_merged$CHR==chr,]
pred_auto <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                         ind.col = df_beta2[["NUM_each_chr"]])

for(chr in 2:22){
  obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
  G   <- obj.bigSNP$genotypes
  
  beta_auto2=beta_auto[dat_merged$CHR==chr,]
  df_beta2=df_beta[dat_merged$CHR==chr,]
  pred_auto2 <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                            ind.col = df_beta2[["NUM_each_chr"]])
  pred_auto=pred_auto+pred_auto2
}

sc <- apply(pred_auto, 2, sd)
# Make TRUE FALSE vector for which converged
convergedAuto <- ((abs(sc - median(sc,na.rm=T)) < 3 * mad(sc,na.rm=T)) & !is.na(sc))
# Average the betas for the models that converged
final_beta_auto <- rowMeans(beta_auto[, convergedAuto])
dat_merged$beta_ldpred_auto=final_beta_auto
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/bipolar/bipolar3/ldpred_results.csv",row.names=F,quote=F)

#####################################perform LDpred - PTSD ###############################

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=fread("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD1/SORTED_PTSD_EA9_AA7_LA1_SA2_ALL_study_specific_PCs1.txt",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("MarkerName","Allele1","Allele2","Effect","StdErr")]
colnames(dat)=c("ID","A1","A2","BETA","SE")
dat$n_eff=9000
colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]


NCORES <- 20
# Open a temporary file
tmp <- tempfile(tmpdir = "/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD1/tmp-data/")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

chr=1
load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
dat_merged2=dat_merged[dat_merged$CHR==chr,]
corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
corr <- as_SFBM(corr0, tmp)
for(chr in 2:22){
  load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
  dat_merged2=dat_merged[dat_merged$CHR==chr,]
  corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
  corr$add_columns(corr0, nrow(corr))
}
save(list="corr",file="/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD1/tmp-data/corr_file_final_PTSD_file1.RData")  

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD1/ldpred_results.csv",row.names=F,quote=F)
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
chr=1
obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
G   <- obj.bigSNP$genotypes

beta_auto2=beta_auto[dat_merged$CHR==chr,]
df_beta2=df_beta[dat_merged$CHR==chr,]
pred_auto <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                         ind.col = df_beta2[["NUM_each_chr"]])

for(chr in 2:22){
  obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
  G   <- obj.bigSNP$genotypes
  
  beta_auto2=beta_auto[dat_merged$CHR==chr,]
  df_beta2=df_beta[dat_merged$CHR==chr,]
  pred_auto2 <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                            ind.col = df_beta2[["NUM_each_chr"]])
  pred_auto=pred_auto+pred_auto2
}

sc <- apply(pred_auto, 2, sd)
# Make TRUE FALSE vector for which converged
convergedAuto <- ((abs(sc - median(sc,na.rm=T)) < 3 * mad(sc,na.rm=T)) & !is.na(sc))
# Average the betas for the models that converged
final_beta_auto <- rowMeans(beta_auto[, convergedAuto])
dat_merged$beta_ldpred_auto=final_beta_auto
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD1/ldpred_results.csv",row.names=F,quote=F)


##############PTSD2#
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=fread("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD2/pts_eur_freeze2_overall.results",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE","Neff")]
dat$OR=log(dat$OR)
colnames(dat)=c("ID","A1","A2","BETA","SE","N_eff")
colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]


NCORES <- 20
# Open a temporary file
tmp <- tempfile(tmpdir = "/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD2/tmp-data/")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

chr=1
load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
dat_merged2=dat_merged[dat_merged$CHR==chr,]
corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
corr <- as_SFBM(corr0, tmp)
for(chr in 2:22){
  load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
  dat_merged2=dat_merged[dat_merged$CHR==chr,]
  corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
  corr$add_columns(corr0, nrow(corr))
}

save(list="corr",file="/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD2/tmp-data/corr_file_final_PTSD_file1.RData")  

#load("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD2/tmp-data/corr_file_final_PTSD_file1.RData")
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD2/ldpred_results.csv",row.names=F,quote=F)
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
chr=1
obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
G   <- obj.bigSNP$genotypes

beta_auto2=beta_auto[dat_merged$CHR==chr,]
df_beta2=df_beta[dat_merged$CHR==chr,]
pred_auto <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                         ind.col = df_beta2[["NUM_each_chr"]])

for(chr in 2:22){
  obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
  G   <- obj.bigSNP$genotypes
  
  beta_auto2=beta_auto[dat_merged$CHR==chr,]
  df_beta2=df_beta[dat_merged$CHR==chr,]
  pred_auto2 <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                            ind.col = df_beta2[["NUM_each_chr"]])
  pred_auto=pred_auto+pred_auto2
}

sc <- apply(pred_auto, 2, sd)
# Make TRUE FALSE vector for which converged
convergedAuto <- ((abs(sc - median(sc,na.rm=T)) < 3 * mad(sc,na.rm=T)) & !is.na(sc))
# Average the betas for the models that converged
final_beta_auto <- rowMeans(beta_auto[, convergedAuto])
dat_merged$beta_ldpred_auto=final_beta_auto
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/PTSD/PTSD2/ldpred_results.csv",row.names=F,quote=F)

####################################perform LDpred - SCZ ###############################

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=read.table("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz1/pgc.scz.full.2012-04.txt",header=T,stringsAsFactors = FALSE,comment.char="!")
dat=dat[,c("snpid","a1","a2","or","se")]
colnames(dat)=c("ID","A1","A2","BETA","SE")
dat$BETA=log(dat$BETA)
dat$n_eff=9000
colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]
h2_est=0.1

NCORES <- 20
# Open a temporary file
tmp <- tempfile(tmpdir = "/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz1/tmp-data/")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

chr=1
load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
dat_merged2=dat_merged[dat_merged$CHR==chr,]
corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
corr <- as_SFBM(corr0, tmp)
for(chr in 2:22){
  load(file=paste("/pl/active/KellerLab/Yongkang/hapmap3/ldfile/tmp-data/corr0_file_chr",chr,".RData",sep=""))  
  dat_merged2=dat_merged[dat_merged$CHR==chr,]
  corr0=corr0[dat_merged2$NUM_each_chr,dat_merged2$NUM_each_chr]
  corr$add_columns(corr0, nrow(corr))
}
save(list="corr",file="/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz1/tmp-data/corr_file_final_scz_file1.RData")  

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz1/ldpred_results.csv",row.names=F,quote=F)
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
chr=1
obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
G   <- obj.bigSNP$genotypes

beta_auto2=beta_auto[dat_merged$CHR==chr,]
df_beta2=df_beta[dat_merged$CHR==chr,]
pred_auto <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                         ind.col = df_beta2[["NUM_each_chr"]])

for(chr in 2:22){
  obj.bigSNP <- snp_attach(paste("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/chr.",chr,".rds",sep=""))
  G   <- obj.bigSNP$genotypes
  
  beta_auto2=beta_auto[dat_merged$CHR==chr,]
  df_beta2=df_beta[dat_merged$CHR==chr,]
  pred_auto2 <- big_prodMat(G, beta_auto2, #ind.row = ind.val,
                            ind.col = df_beta2[["NUM_each_chr"]])
  pred_auto=pred_auto+pred_auto2
}

sc <- apply(pred_auto, 2, sd)
# Make TRUE FALSE vector for which converged
convergedAuto <- ((abs(sc - median(sc,na.rm=T)) < 3 * mad(sc,na.rm=T)) & !is.na(sc))
# Average the betas for the models that converged
final_beta_auto <- rowMeans(beta_auto[, convergedAuto])
dat_merged$beta_ldpred_auto=final_beta_auto
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz1/ldpred_results.csv",row.names=F,quote=F)

#####SCZ 2

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=fread("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz2/daner_PGC_SCZ52_0513a.hq2",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE","ngt")]
colnames(dat)=c("ID","A1","A2","BETA","SE","n_eff")
dat$BETA=log(dat$BETA)
colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz2/summary_file.csv",row.names=F,col.names=T,quote=F)
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]
h2_est=0.1

NCORES <- 20


#####SCZ 3

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=fread("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz3/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE")]
colnames(dat)=c("ID","A1","A2","BETA","SE")
dat$BETA=log(dat$BETA)
dat$n_eff=30000

colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz3/summary_file.csv",row.names=F,col.names=T,quote=F)
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]
h2_est=0.1

NCORES <- 20


#####SCZ 4

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=fread("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz4/daner_natgen_eas_eur_autosomes",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE","Neff")]
colnames(dat)=c("ID","A1","A2","BETA","SE","n_eff")
dat$BETA=log(dat$BETA)

colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/scz/scz4/summary_file.csv",row.names=F,col.names=T,quote=F)
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]
h2_est=0.1

NCORES <- 20

############Major depression ##############


library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

#dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,stringsAsFactors = FALSE)
#dat=as.data.frame(dat)
#dat2=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,fill=TRUE,stringsAsFactors = FALSE)
#dat2=as.data.frame(dat2)
dat3=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,stringsAsFactors = FALSE,nrows=9533408)
dat3=as.data.frame(dat3)
dat4=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=F,stringsAsFactors = FALSE,skip=9533409,nrow=9533470-9533409)
dat4=as.data.frame(dat4)
dat5=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=F,stringsAsFactors = FALSE,skip=9533470)
dat5=as.data.frame(dat5)
colnames(dat5)=colnames(dat3)
dat=rbind(dat3,dat5)

dat=dat[,c("SNP","A1","A2","OR","SE","Neff")]
colnames(dat)=c("ID","A1","A2","BETA","SE","n_eff")
dat$BETA=log(dat$BETA)

colnames(dat)[1]="rsid"
hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
dir.create("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/major_depression")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/major_depression/summary_file.csv",row.names=F,col.names=T,quote=F)
df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]
ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]

#######OCD ########
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_ocd/ocd_aug2017",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE")]
colnames(dat)=c("ID","A1","A2","BETA","SE")
dat$n_eff=30000
dat$BETA=log(dat$BETA)

colnames(dat)[1]="rsid"

hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
dir.create("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/ocd")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/ocd/summary_file.csv",row.names=F,col.names=T,quote=F)

######AUTISM######
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

#dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,stringsAsFactors = FALSE)
#dat=as.data.frame(dat)
#dat2=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,fill=TRUE,stringsAsFactors = FALSE)
#dat2=as.data.frame(dat2)
dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_autism/iPSYCH-PGC_ASD_Nov2017",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE")]
colnames(dat)=c("ID","A1","A2","BETA","SE")
dat$n_eff=30000
dat$BETA=log(dat$BETA)

colnames(dat)[1]="rsid"

hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
dir.create("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/autism")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/autism/summary_file.csv",row.names=F,col.names=T,quote=F)


######anxiety#######

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

#dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,stringsAsFactors = FALSE)
#dat=as.data.frame(dat)
#dat2=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,fill=TRUE,stringsAsFactors = FALSE)
#dat2=as.data.frame(dat2)
dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_anxiety/anxiety.meta.full.fs.tbl",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNPID","Allele1","Allele2","Effect","StdErr","TotalN")]
colnames(dat)=c("ID","A1","A2","BETA","SE","n_eff")
#dat$n_eff=30000
#dat$BETA=log(dat$BETA)

colnames(dat)[1]="rsid"

hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
dir.create("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/anxiety_fs")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/anxiety_fs/summary_file.csv",row.names=F,col.names=T,quote=F)


dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_anxiety/anxiety.meta.full.cc.tbl",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNPID","Allele1","Allele2","Effect","StdErr","TotalN")]
colnames(dat)=c("ID","A1","A2","BETA","SE","n_eff")
#dat$n_eff=30000
#dat$BETA=log(dat$BETA)

colnames(dat)[1]="rsid"

hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
dir.create("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/anxiety_cc")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/anxiety_cc/summary_file.csv",row.names=F,col.names=T,quote=F)

######ADHD#######

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

#dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,stringsAsFactors = FALSE)
#dat=as.data.frame(dat)
#dat2=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,fill=TRUE,stringsAsFactors = FALSE)
#dat2=as.data.frame(dat2)
dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/PGC_ADHD/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE","Neff")]
colnames(dat)=c("ID","A1","A2","BETA","SE","n_eff")
#dat$n_eff=30000
dat$BETA=log(dat$BETA)

colnames(dat)[1]="rsid"

hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
dir.create("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/ADHD")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/ADHD/summary_file.csv",row.names=F,col.names=T,quote=F)



######Neuroticism#### Don't use it



library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

#dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,stringsAsFactors = FALSE)
#dat=as.data.frame(dat)
#dat2=fread("/pl/active/KellerLab/Yongkang/PGS_web/pgc_major_depression/PGC_non_UKB_depression.txt",header=T,fill=TRUE,stringsAsFactors = FALSE)
#dat2=as.data.frame(dat2)
dat=fread("/pl/active/KellerLab/Yongkang/PGS_web/neuroticism_Mats/sumstats_neuroticism_ctg_format.txt",header=T,stringsAsFactors = FALSE)
dat=as.data.frame(dat)
dat=dat[,c("SNP","A1","A2","OR","SE","Neff")]
colnames(dat)=c("ID","A1","A2","BETA","SE","n_eff")
#dat$n_eff=30000
dat$BETA=log(dat$BETA)

colnames(dat)[1]="rsid"

hapmap_info=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
hapmap_info=as.data.frame(hapmap_info)

hapmap_info2=hapmap_info[which(!is.na(hapmap_info$MENG_ID)),]
hapmap_info2=hapmap_info2[which(hapmap_info2$MENG_ID!=""),]
hapmap_info2=hapmap_info2[,c("rsid","MENG_ID","chr","GRCH38pos")]
dat_merged=merge(hapmap_info2,dat,by="rsid")

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile2=mapfile[match(dat_merged$MENG_ID,mapfile$ID),]
dat_merged$BETA[mapfile2$A1!=dat_merged$A1]=-dat_merged$BETA[mapfile2$A1!=dat_merged$A1]
dat_merged=cbind(dat_merged[,c(2,7,8,9)],mapfile2[,c(1,3,4,5,6,7,8)])
dat_merged=dat_merged[,c(1,5,6,7,8,2,3,4,9,10,11)]
dat_merged=dat_merged[order(dat_merged$"_NUM_ID"),]
colnames(dat_merged)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
dir.create("/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/neuroticism")
fwrite(dat_merged,"/pl/active/KellerLab/Yongkang/new_ldpred_results/pgc/neuroticism/summary_file.csv",row.names=F,col.names=T,quote=F)

NCORES <- 20
######################Generate new pgs info######
chr=1
map=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/new_LD/indep_UKB/ldfile/tmp-data/mapfile_chr",chr,".txt",sep=""),header=T,stringsAsFactors=FALSE)
map=as.data.frame(map)
aa=seq(1,nrow(map),length.out=51)
aa=as.integer(aa)
aa[1]=0
dir.create("/pl/active/KellerLab/Yongkang/hapmap3/split_chr")
for(iter in 1:50){
  map2=map[(aa[iter]+1):aa[iter+1],]
  fwrite(map2,paste("/pl/active/KellerLab/Yongkang/hapmap3/split_chr/mapfile_chr",chr,".idx",iter,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  write.table(map2$rsid,paste("/pl/active/KellerLab/Yongkang/hapmap3/split_chr/snplist_chr",chr,".idx",iter,".txt",sep=""),row.names=F,col.names=F,quote=F)
}

for(chr in 2:22){
  map=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/new_LD/indep_UKB/ldfile/tmp-data/mapfile_chr",chr,".txt",sep=""),header=T,stringsAsFactors=FALSE)
  map=as.data.frame(map)
  aa=seq(1,nrow(map),length.out=51)
  aa=as.integer(aa)
  aa[1]=0
  for(iter in 1:50){
    map2=map[(aa[iter]+1):aa[iter+1],]
    fwrite(map2,paste("/pl/active/KellerLab/Yongkang/hapmap3/split_chr/mapfile_chr",chr,".idx",iter,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
    write.table(map2$rsid,paste("/pl/active/KellerLab/Yongkang/hapmap3/split_chr/snplist_chr",chr,".idx",iter,".txt",sep=""),row.names=F,col.names=F,quote=F)
  }
  
}

###Divide each chromosome data to 50 subsets for reducing the time consuming.
ml gcc
chr=1;idx=1
/projects/yoki5348/vcftools_0.1.13/bin/vcftools --gzvcf /pl/active/KellerLab/Yongkang/hapmap3/indep_samp_hap/UKB_chr"$chr".indep.vcf.recode.vcf --snps /pl/active/KellerLab/Yongkang/hapmap3/split_chr/snplist_chr"$chr".idx"$idx".txt --recode --recode-INFO-all --out /scratch/summit/yoki5348/indep_samp_hap/UKB_chr"$chr".idx"$idx".indep.vcf
chr=1
for idx in {1..50}
do
/projects/yoki5348/vcftools_0.1.13/bin/vcftools --gzvcf /pl/active/KellerLab/Yongkang/hapmap3/indep_samp_hap/UKB_chr"$chr".indep.vcf.recode.vcf --snps /pl/active/KellerLab/Yongkang/hapmap3/split_chr/snplist_chr"$chr".idx"$idx".txt --recode --recode-INFO-all --out /scratch/summit/yoki5348/indep_samp_hap/UKB_chr"$chr".idx"$idx".indep.vcf
done

dat=fread("/pl/active/KellerLab/Yongkang/hapmap3/hapmap3_liftover_results.txt",header=T,stringsAsFactors=FALSE)
dat=as.data.frame(dat)


