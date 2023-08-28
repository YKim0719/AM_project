#singularity run /projects/lessem/singularity/openmx.sif Rscript /pl/active/KellerLab/Yongkang/AM_source_code/SEM_for_AM_with_internal_summary_statistics_WHR_w.o.BMI_w.o.invNorm_each_sex.R
####Check regression results #######3
library(data.table)

source("/pl/active/KellerLab/Yongkang/AM_source_code/source_code_for_SEM.R")
source("/pl/active/KellerLab/Yongkang/AM_source_code/source_code_for_various_hpgs_handle_ver2.R")
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno2=fread("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/pheno_w.o.covariate_adjustment_european_ver2.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)
pheno2$IID=pheno2$FID
pheno3=fread("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/pheno_w.o.covariate_adjustment_test.txt",header=T,stringsAsFactors=FALSE)
pheno3=as.data.frame(pheno3)
pheno3$IID=pheno3$FID
pheno2=pheno2[match(pheno3$IID,pheno2$IID),]

binary=c()
for(idx in 3:ncol(pheno2)){
  if(length(table(pheno2[,idx]))<=2){
    binary=c(binary,idx)
  }
}

conti=c(3:15,18,19,22,23,26,28,29,33:35)
pheno_invnorm=fread("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/pheno_w.invNorm_first_european_test_w.o.adjustment.txt",header=T,stringsAsFactors=FALSE)
pheno_invnorm=as.data.frame(pheno_invnorm)
pheno_invnorm$IID=pheno_invnorm$FID
pheno_invnorm=pheno_invnorm[match(pheno2$IID,pheno_invnorm$IID),]
marriage_year=fread("/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/Assortment_Project/Marriage_Years/UKB_SpousePairs_MarriageYears.txt",header=T,stringsAsFactors=FALSE)
marriage_year=as.data.frame(marriage_year)
marriage_year=marriage_year[match(pheno2$IID,marriage_year$f.eid),]
marriage_year=marriage_year[,c(1,3,2,4)]

covar=fread("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/covar_file_european_BMI.txt",header=T,stringsAsFactors=FALSE)
covar=as.data.frame(covar)
covar=covar[match(pheno2$IID,covar$IID),]
covar_EA=fread("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/covar_file_european_income.txt",header=T,stringsAsFactors=FALSE)
covar_EA=as.data.frame(covar_EA)
covar_EA=covar[match(pheno2$IID,covar_EA$IID),]
pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
pheno_invnorm=pheno_invnorm[match(pheno2$IID,pheno_invnorm$IID),]
covar=covar[which(marriage_year$Participant.Marriage.Year<=1999),]
covar_EA=covar_EA[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]

pheno_invnorm_after4=fread("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/pheno_w.invNorm_later_european_test_w.center_w.PC.txt",header=T,stringsAsFactors=FALSE)
pheno_invnorm_after4=as.data.frame(pheno_invnorm_after4)
pheno_invnorm_after4=pheno_invnorm_after4[match(pheno2$IID,pheno_invnorm_after4$IID),]


spo=fread("/pl/active/KellerLab/Yongkang/UKB_spo_pair.csv",header=T,stringsAsFactors=FALSE)
spo=as.data.frame(spo)

for(idx2 in conti){
  phenoname=colnames(pheno2)[idx2]
  
  pheno_invnorm_after4_test=pheno_invnorm_after4[,c("IID",phenoname)]
  pheno_invnorm="pheno_invnorm_after4"
  
  dir.create(paste("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/regenie/SEM_results_no_trend/",sep="") )
  dir_name=paste("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/regenie/SEM_results_no_trend/",phenoname,sep="") 
  time1=Sys.time()
  dat=fread(paste("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/regenie/",phenoname,"/w.o.invNorm_w.PC_w.center/individual_hpgs_merged.csv",sep=""),header=T,stringsAsFactors=FALSE)
  dat=as.data.frame(dat)
  dat=dat[match(pheno2$IID,dat$IID),]
  k=2*((var(dat[,2],na.rm=TRUE)+var(dat[,3],na.rm=TRUE))/2-(cov(dat[,2],dat[,3],use="complete.obs")))
  dat[,2:3]=dat[,2:3]/sqrt(k)
  
  dat_fat=dat[match(spo$IID_Fat,dat$IID),]
  dat_mot=dat[match(spo$IID_Mot,dat$IID),]
  
  exclude=which(is.na(dat[,2])|is.na(get(paste(pheno_invnorm,"_test",sep=""))[,2]))
  if(length(exclude)!=0){
    dat2=dat[-exclude,]
    marriage_year=marriage_year[-exclude,]
    pheno_used=get(paste(pheno_invnorm,"_test",sep=""))[-exclude,]
  }else{
    dat2=dat
    pheno_used=get(paste(pheno_invnorm,"_test",sep=""))
  }
  pheno_used=cbind(pheno_used$IID,pheno_used)
  colnames(pheno_used)=c("FID","IID",phenoname)
  pheno_fat=pheno_used[match(spo$IID_Fat,dat2$IID),]
  pheno_mot=pheno_used[match(spo$IID_Mot,dat2$IID),]
  dir.create(dir_name)
  spo_fat=NA;spo_mot=NA;bin_idx=NA
  
  
  pgs=dat2[,2]+dat2[,3]
  lmfit=lm(pheno_used[,3]~pgs)
  R2_ref=summary(lmfit)$r.squared
  g_result=try(perform_SEM_gc_gt_differ(phenoname,dat2,pheno_used,R2_ref,spo$IID_Fat,spo$IID_Mot,dir_name))
  global_g_result=try(perform_SEM_gc_gt_same(phenoname,dat2,pheno_used,R2_ref,spo$IID_Fat,spo$IID_Mot,dir_name))
  cor_ref=cor(pheno_fat[,2],pheno_mot[,2],use="complete.obs")
  
  bb=mxCompare(g_result,global_g_result )
  write.csv(bb,paste(dir_name,"/lrt_test_results.csv",sep=""),row.names=F,quote=F)
  
  
  time2=Sys.time()
  print(time2-time1)
}

binary=binary[-4]
for(idx2 in binary){
  phenoname=colnames(pheno2)[idx2]
  pheno_binary=pheno2[,c("IID",phenoname)]  
  pheno_binary[,2]=pheno_binary[,2]-1
  if(min(table(pheno_binary[,2]))>=5000){
  dir.create(paste("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/regenie/SEM_results_no_trend/",sep="") )
  dir_name=paste("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/regenie/SEM_results_no_trend/",phenoname,sep="") 
  dat=fread(paste("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/regenie/",phenoname,"/w.o.invNorm_w.PC_w.center/individual_hpgs_merged.csv",sep=""),header=T,stringsAsFactors=FALSE)
  dat=as.data.frame(dat)
  dat=dat[match(pheno2$IID,dat$IID),]
  k=2*((var(dat[,2],na.rm=TRUE)+var(dat[,3],na.rm=TRUE))/2-(cov(dat[,2],dat[,3],use="complete.obs")))
  dat[,2:3]=dat[,2:3]/sqrt(k)
  
  dat_fat=dat[match(spo$IID_Fat,dat$IID),]
  dat_mot=dat[match(spo$IID_Mot,dat$IID),]
  
  exclude=which(is.na(dat[,2])|is.na(pheno_binary[,2]))
  if(length(exclude)!=0){
    dat2=dat[-exclude,]
    marriage_year=marriage_year[-exclude,]
    pheno_used=pheno_binary[-exclude,]
  }else{
    dat2=dat
    pheno_used=pheno_binary
  }
  pheno_used=cbind(pheno_used$IID,pheno_used)
  colnames(pheno_used)=c("FID","IID",phenoname)
  pheno_fat=pheno_used[match(spo$IID_Fat,dat2$IID),]
  pheno_mot=pheno_used[match(spo$IID_Mot,dat2$IID),]
  dir.create(dir_name)
  cor_ref=read.csv("/pl/active/KellerLab/Yongkang/AM_source_code/UKB_pheno_r2_cor_ver3.csv",header=T,stringsAsFactors=FALSE)
  if(phenoname=="SmkInit"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="SmokingInitiation_gscan"]}
  if(phenoname=="SmkCes"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="SmokingCessation_gscan"]}
  if(phenoname=="icd_mdd"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="pgs_major_depression_non_UKB"]}
  if(phenoname=="life_mdd"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="pgs_major_depression_non_UKB"]}
  if(phenoname=="bdd"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="pgs_major_depression_non_UKB"]}
  if(phenoname=="nervous"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="anxiety_cc"]}
  if(phenoname%in%c("snore","fed_up","T2D")){
    pheno_val=pheno2[,phenoname]
    pheno_val=pheno_val-1
    lmfit=glm(pheno_val~.,data=covar,family=binomial)
    fitted=predict(lmfit,newdata=covar,type="response")
    pheno_val=scale(pheno_val-fitted)
    pheno_dat=data.frame(IID=pheno2$IID,pheno_val=pheno_val)
    colnames(pheno_dat)[2]=phenoname
    pgs=dat[,2]+dat[,3]
    lmfit=lm(pheno_dat[,2]~pgs)
    R2_ref=summary(lmfit)$r.squared
  }

     g_result=try(perform_SEM_gc_gt_differ_binary(phenoname,dat2,pheno_used,R2_ref,spo$IID_Fat,spo$IID_Mot,dir_name))
    global_g_result=try(perform_SEM_gc_gt_same_binary(phenoname,dat2,pheno_used,R2_ref,spo$IID_Fat,spo$IID_Mot,dir_name))
  bb=mxCompare(g_result,global_g_result )
  write.csv(bb,paste(dir_name,"/lrt_test_results.csv",sep=""),row.names=F,quote=F)
  
}else{
  dir.create(paste("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/regenie/SEM_results_no_trend/",sep="") )
  dir_name=paste("/pl/active/KellerLab/Yongkang/UKB_internal_with_various_adjustment/regenie/SEM_results_no_trend/",phenoname,sep="") 
  dat=read.csv(paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/hpgs_collection/",phenoname,"/hpgs_pc25.csv",sep=""),header=T,stringsAsFactors=FALSE)
  dat=dat[match(pheno2$IID,dat$IID),]
  k=2*((var(dat[,2],na.rm=TRUE)+var(dat[,3],na.rm=TRUE))/2-(cov(dat[,2],dat[,3],use="complete.obs")))
  dat[,2:3]=dat[,2:3]/sqrt(k)
  
  dat=dat[match(pheno2$IID,dat$IID),]
  k=2*((var(dat[,2],na.rm=TRUE)+var(dat[,3],na.rm=TRUE))/2-(cov(dat[,2],dat[,3],use="complete.obs")))
  dat[,2:3]=dat[,2:3]/sqrt(k)
  
  dat_fat=dat[match(spo$IID_Fat,dat$IID),]
  dat_mot=dat[match(spo$IID_Mot,dat$IID),]
  
  exclude=which(is.na(dat[,2])|is.na(pheno_binary[,2]))
  if(length(exclude)!=0){
    dat2=dat[-exclude,]
    marriage_year=marriage_year[-exclude,]
    pheno_used=pheno_binary[-exclude,]
  }else{
    dat2=dat
    pheno_used=pheno_binary
  }
  pheno_used=cbind(pheno_used$IID,pheno_used)
  colnames(pheno_used)=c("FID","IID",phenoname)
  pheno_fat=pheno_used[match(spo$IID_Fat,dat2$IID),]
  pheno_mot=pheno_used[match(spo$IID_Mot,dat2$IID),]
  dir.create(dir_name)
  cor_ref=read.csv("/pl/active/KellerLab/Yongkang/AM_source_code/UKB_pheno_r2_cor_ver3.csv",header=T,stringsAsFactors=FALSE)
  if(phenoname=="SmkInit"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="SmokingInitiation_gscan"]}
  if(phenoname=="SmkCes"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="SmokingCessation_gscan"]}
  if(phenoname=="icd_mdd"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="pgs_major_depression_non_UKB"]}
  if(phenoname=="life_mdd"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="pgs_major_depression_non_UKB"]}
  if(phenoname=="bdd"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="pgs_major_depression_non_UKB"]}
  if(phenoname=="nervous"){ R2_ref=cor_ref$R2_ref[cor_ref$pheno=="anxiety_cc"]}
  if(phenoname%in%c("snore","fed_up","T2D")){
    pheno_val=pheno2[,phenoname]
    pheno_val=pheno_val-1
    lmfit=glm(pheno_val~.,data=covar,family=binomial)
    fitted=predict(lmfit,newdata=covar,type="response")
    pheno_val=scale(pheno_val-fitted)
    pheno_dat=data.frame(IID=pheno2$IID,pheno_val=pheno_val)
    colnames(pheno_dat)[2]=phenoname
    pgs=dat[,2]+dat[,3]
    lmfit=lm(pheno_dat[,2]~pgs)
    R2_ref=summary(lmfit)$r.squared
  }
  
  g_result=try(perform_SEM_gc_gt_differ_w.o.Y(phenoname,dat2,R2_ref,spo$IID_Fat,spo$IID_Mot,dir_name))
  global_g_result=try(perform_SEM_gc_gt_same_w.o.Y(phenoname,dat2,R2_ref,spo$IID_Fat,spo$IID_Mot,dir_name))
  

  bb=mxCompare(g_result,global_g_result )
  write.csv(bb,paste(dir_name,"/lrt_test_results.csv",sep=""),row.names=F,quote=F)
}
}

