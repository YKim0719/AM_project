library(data.table)
cor_inf=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_r2_cor_ver2.csv",header=T,stringsAsFactors=FALSE)
cor_inf=as.data.frame(cor_inf)
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_test.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_phenofiles_no_scaled_addmore.csv",header=T,stringsAsFactors=FALSE)

marriage_year=fread("/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/Assortment_Project/Marriage_Years/UKB_SpousePairs_MarriageYears.txt",header=T,stringsAsFactors=FALSE)
marriage_year=as.data.frame(marriage_year)
marriage_year=marriage_year[match(pheno2$IID,marriage_year$f.eid),]
marriage_year=marriage_year[,c(1,3,2,4)]

covar=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_covarfile3.txt",header=T,stringsAsFactors=FALSE)
covar=as.data.frame(covar)
covar=covar[match(pheno2$IID,covar$id),]
covar=covar[,-c(4:13)]

spo=fread("/pl/active/KellerLab/Yongkang/UKB_spo_pair.csv",header=T,stringsAsFactors=FALSE)
spo=as.data.frame(spo)

bin_range_year=3
bin_range_std=5
bin_year=TRUE; bin_std=TRUE

eigen_vec=fread("/pl/active/KellerLab/Yongkang/UKB_PC_distributed.csv",header=T,stringsAsFactors=FALSE)
eigen_vec=as.data.frame(eigen_vec)
eigen_vec=eigen_vec[,1:41]
no.pc=40
dir.create("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_internal/SEM_analysis_results/")
eigen_vec2=eigen_vec[match(pheno2$IID,eigen_vec[,1]),1:(no.pc+1)]

covar2=covar[,-1]
covar2=cbind(covar2,eigen_vec2[,-1])
covar2$year_sex=covar2$Year_of_birth*covar2$Sex

pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
covar2=covar2[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]

phenonames=colnames(pheno2)
idx=6
phenoname=phenonames[idx]
if(length(unique(na.omit(pheno2[,phenoname])))==2){
  pheno_val=pheno2[,phenoname]
  lmfit=glm(pheno_val~.,data=covar2,family=binomial)
  fitted=predict(lmfit,newdata=covar2,type="response")
  pheno_val=scale(pheno_val-fitted)
  pheno_dat=data.frame(IID=pheno2$IID,pheno_val=pheno_val)
  colnames(pheno_dat)[2]=phenoname
}else{
  pheno_val=pheno2[,phenoname]
  lmfit=lm(pheno_val~.,data=covar2)
  fitted=predict(lmfit,newdata=covar2)
  pheno_val=scale(pheno_val-fitted)
  pheno_dat=data.frame(IID=pheno2$IID,pheno_val=pheno_val)
  colnames(pheno_dat)[2]=phenoname
  
}

hpgs=read.csv(paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_internal/hpgs_collection/",phenoname,"/hpgs_pc25.csv",sep=""),header=T,stringsAsFactors=FALSE)

hpgs2=hpgs[match(pheno2$IID,hpgs$IID),]

exclude=which(is.na(hpgs2$pgs_inf1)|is.na(pheno2[,phenoname])|is.na(pheno_dat[,2]))
if(length(exclude)!=0){
  hpgs2=hpgs2[-exclude,]
  covar2=covar2[-exclude,]
  marriage_year=marriage_year[-exclude,]
  pheno2=pheno2[-exclude,]
  pheno_dat=pheno_dat[-exclude,]
}
pheno_fat=pheno_dat[match(spo$IID_Fat,pheno_dat[,1]),]
pheno_mot=pheno_dat[match(spo$IID_Mot,pheno_dat[,1]),]

pgs=hpgs2[,2]+hpgs2[,3]
lmfit=lm(pheno_dat[,2]~pgs)
results=data.frame(phenoname=phenoname,Rsq=summary(lmfit)$r.squared,cor_UKB=cor(pheno_fat[,2],pheno_mot[,2],use="complet"))

for(idx in 7:ncol(pheno2)){
  phenoname=phenonames[idx]
  if(length(unique(na.omit(pheno2[,phenoname])))==2){
    pheno_val=pheno2[,phenoname]
    lmfit=glm(pheno_val~.,data=covar2,family=binomial)
    fitted=predict(lmfit,newdata=covar2,type="response")
    pheno_val=scale(pheno_val-fitted)
    pheno_dat=data.frame(IID=pheno2$IID,pheno_val=pheno_val)
    colnames(pheno_dat)[2]=phenoname
  }else{
    pheno_val=pheno2[,phenoname]
    lmfit=lm(pheno_val~.,data=covar2)
    fitted=predict(lmfit,newdata=covar2)
    pheno_val=scale(pheno_val-fitted)
    pheno_dat=data.frame(IID=pheno2$IID,pheno_val=pheno_val)
    colnames(pheno_dat)[2]=phenoname
    
  }
  
  hpgs=try(read.csv(paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_internal/hpgs_collection/",phenoname,"/hpgs_pc25.csv",sep=""),header=T,stringsAsFactors=FALSE))
  if(class(hpgs)[1]!="try-error"){
    hpgs2=hpgs[match(pheno2$IID,hpgs$IID),]
    
    exclude=which(is.na(hpgs2$pgs_inf1)|is.na(pheno2[,phenoname])|is.na(pheno_dat[,2]))
    if(length(exclude)!=0){
      hpgs2=hpgs2[-exclude,]
      covar2=covar2[-exclude,]
      marriage_year=marriage_year[-exclude,]
      pheno2=pheno2[-exclude,]
      pheno_dat=pheno_dat[-exclude,]
    }
    
    pheno_fat=pheno_dat[match(spo$IID_Fat,pheno_dat[,1]),]
    pheno_mot=pheno_dat[match(spo$IID_Mot,pheno_dat[,1]),]
    
    pgs=hpgs2[,2]+hpgs2[,3]
    lmfit=lm(pheno_dat[,2]~pgs)
    results2=data.frame(phenoname=phenoname,Rsq=summary(lmfit)$r.squared,cor_UKB=cor(pheno_fat[,2],pheno_mot[,2],use="complet"))
    results=rbind(results,results2)
  }  
 
}

cor_inf2=cor_inf[match(results[,1],cor_inf[,1]),]
cor_inf2[,1]=results[,1]
cor_inf2$R2_UKB=results[,2]
cor_inf2$cor_UKB=results[,3]
cor_inf2=cor_inf2[,-c(6,7)]
write.csv(cor_inf2,"/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_internal/R_sqaured_internal.csv",row.names=F,quote=F)



mkdir /pl/active/KellerLab/Yongkang/UK_internal
cp -r /rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_internal /pl/active/KellerLab/Yongkang/UK_internal
