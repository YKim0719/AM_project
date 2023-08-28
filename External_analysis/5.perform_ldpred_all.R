args <- commandArgs(trailingOnly =TRUE)
idx=args[1]
dir_name_input=args[1]
dir_name_output=args[2]

#Total 45 names here
#dir.create("/pl/active/KellerLab/Yongkang/ldsc_results/")
#dirs=dir("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/")
#aa=data.frame(input=paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",dirs,sep=""),output=paste("/pl/active/KellerLab/Yongkang/ldsc_results/",dirs,sep=""))
#write.csv(aa,"/pl/active/KellerLab/Yongkang/ldsc_results/id_collection.csv",row.names=F,quote=F)

aa=read.csv("/pl/active/KellerLab/Yongkang/ldsc_results/id_collection.csv",header=T,stringsAsFactors=FALSE)
dir_name_input=aa[idx,1]
dir_name_output=aa[idx,2]

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

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


dat=fread(paste(dir_name_input,"/sumstat_pheno_all.txt",sep=""),header=T,stringsAsFactors=FALSE)
dat=as.data.frame(dat)
dat$n_eff=500000
mapfile2=mapfile[match(dat$snpid,mapfile$ID),]
used_var=which((mapfile2$A0==dat$a1&mapfile2$A1==dat$a2)|(mapfile2$A0==dat$a2&mapfile2$A1==dat$a1))
mapfile2=mapfile2[used_var,]
dat_merged=dat[used_var,]
dat_merged$beta[mapfile2$A1!=dat_merged$a2]=-dat_merged$beta[mapfile2$A1!=dat_merged$a2]
dat_merged=cbind(dat_merged[,c("snpid","beta","se","n_eff")],mapfile2[,c(1,3,4,5,6,7,8)])
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
if(h2_est<0 | h2_est>1){h2_est=0.1}


NCORES <- 10
# Open a temporary file
dir.create(dir_name_output)
tmp <- tempfile(tmpdir = paste(dir_name_output,"/tmp-data/",sep=""))
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
save(list="corr",file=paste(dir_name_output,"/tmp-data/corr_file_final.RData",sep=""))  

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
fwrite(dat_merged,paste(dir_name_output,"/ldpred_results.csv",sep=""),row.names=F,quote=F)
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
fwrite(dat_merged,paste(dir_name_output,"/ldpred_results.csv",sep=""),row.names=F,quote=F)

system(paste("rm -r ",dir_name_output,"/tmp-data"),sep="")