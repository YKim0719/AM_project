args <- commandArgs(trailingOnly =TRUE)
idx=as.numeric(args[1])
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)
library(stringr)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british_include_new_alz.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)
phenoname=colnames(pheno2)[idx]

mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile=as.data.frame(mapfile)
if(phenoname=="alz"){
  #paste("/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/alz/chr5_train_alz.regenie",sep="")
  chr=1
  dat=fread(paste("/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/alz/chr",chr,"_train_",phenoname,".regenie",sep=""),header=T,stringsAsFactors=FALSE)
  dat=as.data.frame(dat)
  mapfile2=mapfile[match(dat$ID,mapfile$ID),]
  
  check=which((mapfile2$A1==dat$ALLELE1) &(mapfile2$A0==dat$ALLELE0))
  check2=which((mapfile2$A0==dat$ALLELE1) &(mapfile2$A1==dat$ALLELE0))
  check=c(check,check2)
  dat=dat[check,]
  mapfile2=mapfile[match(dat$ID,mapfile$ID),]
  
  dat$BETA[mapfile2$A1!=dat$ALLELE1]=-dat$BETA[mapfile2$A1!=dat$ALLELE1]
  dat=cbind(dat[,c("ID","CHROM","GENPOS")],mapfile2[,c("A0","A1")],dat[,c("BETA","SE","N")],mapfile2[,c("ld","_NUM_ID","NUM_each_chr")])
  dat=dat[order(dat$"_NUM_ID"),]
  colnames(dat)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
  dat_merged=dat
  for(chr in 2:22){
    dat=fread(paste("/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/alz/chr",chr,"_train_",phenoname,".regenie",sep=""),header=T,stringsAsFactors=FALSE)
    dat=as.data.frame(dat)
    mapfile2=mapfile[match(dat$ID,mapfile$ID),]
    check=which((mapfile2$A1==dat$ALLELE1) &(mapfile2$A0==dat$ALLELE0))
    check2=which((mapfile2$A0==dat$ALLELE1) &(mapfile2$A1==dat$ALLELE0))
    check=c(check,check2)
    dat=dat[check,]
    mapfile2=mapfile[match(dat$ID,mapfile$ID),]
    
    dat$BETA[mapfile2$A1!=dat$ALLELE1]=-dat$BETA[mapfile2$A1!=dat$ALLELE1]
    dat=cbind(dat[,c("ID","CHROM","GENPOS")],mapfile2[,c("A0","A1")],dat[,c("BETA","SE","N")],mapfile2[,c("ld","_NUM_ID","NUM_each_chr")])
    dat=dat[order(dat$"_NUM_ID"),]
    colnames(dat)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
    dat_merged=rbind(dat_merged,dat)
    
  }
}else{
  chr=1
  dat=fread(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"/train_chr",chr,"_",phenoname,".regenie",sep=""),header=T,stringsAsFactors=FALSE)
  dat=as.data.frame(dat)
  mapfile2=mapfile[match(dat$ID,mapfile$ID),]
  
  check=which((mapfile2$A1==dat$ALLELE1) &(mapfile2$A0==dat$ALLELE0))
  check2=which((mapfile2$A0==dat$ALLELE1) &(mapfile2$A1==dat$ALLELE0))
  check=c(check,check2)
  dat=dat[check,]
  mapfile2=mapfile[match(dat$ID,mapfile$ID),]
  
  dat$BETA[mapfile2$A1!=dat$ALLELE1]=-dat$BETA[mapfile2$A1!=dat$ALLELE1]
  dat=cbind(dat[,c("ID","CHROM","GENPOS")],mapfile2[,c("A0","A1")],dat[,c("BETA","SE","N")],mapfile2[,c("ld","_NUM_ID","NUM_each_chr")])
  dat=dat[order(dat$"_NUM_ID"),]
  colnames(dat)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
  dat_merged=dat
  for(chr in 2:22){
    dat=fread(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"/train_chr",chr,"_",phenoname,".regenie",sep=""),header=T,stringsAsFactors=FALSE)
    dat=as.data.frame(dat)
    mapfile2=mapfile[match(dat$ID,mapfile$ID),]
    check=which((mapfile2$A1==dat$ALLELE1) &(mapfile2$A0==dat$ALLELE0))
    check2=which((mapfile2$A0==dat$ALLELE1) &(mapfile2$A1==dat$ALLELE0))
    check=c(check,check2)
    dat=dat[check,]
    mapfile2=mapfile[match(dat$ID,mapfile$ID),]
    
    dat$BETA[mapfile2$A1!=dat$ALLELE1]=-dat$BETA[mapfile2$A1!=dat$ALLELE1]
    dat=cbind(dat[,c("ID","CHROM","GENPOS")],mapfile2[,c("A0","A1")],dat[,c("BETA","SE","N")],mapfile2[,c("ld","_NUM_ID","NUM_each_chr")])
    dat=dat[order(dat$"_NUM_ID"),]
    colnames(dat)=c("ID","CHR","POS","A0","A1","beta","beta_se","n_eff","ld","_NUM_ID_","NUM_each_chr")
    dat_merged=rbind(dat_merged,dat)
    
  }
  
}


df_beta=dat_merged[,c("beta", "beta_se", "n_eff", "_NUM_ID_","NUM_each_chr")]

ld=dat_merged[,"ld"]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]
if(h2_est<0){
  h2_est=0.001
}

NCORES <- 20
# Open a temporary file
dir.create(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,sep=""))
if(sum(dir(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,sep=""))=="corr_file_final_bipolar_file1.RData")!=1){
  tmp <- tempfile(tmpdir = paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,sep=""))
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
  save(list="corr",file=paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"/corr_file_final_bipolar_file1.RData",sep=""))  
  
}else{
  load(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"/corr_file_final_bipolar_file1.RData",sep=""))
}

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
fwrite(dat_merged,paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"/ldpred_results.csv",sep=""),row.names=F,quote=F)

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
fwrite(dat_merged,paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"/ldpred_results.csv",sep=""),row.names=F,quote=F)

system(paste("rm ",tmp, ".sbk",sep=""))
system(paste("rm /pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"/corr_file_final_bipolar_file1.RData",sep=""))
