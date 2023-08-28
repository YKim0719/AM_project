args <- commandArgs(trailingOnly =TRUE)
idx=as.numeric(args[1])

library(data.table)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british_include_new_alz.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)
phenonames=dir("/pl/active/KellerLab/Yongkang/ldsc_results/")
cor_inf=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_r2_cor_ver2.csv",header=T,stringsAsFactors=FALSE)
cor_inf=as.data.frame(cor_inf)

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
eigen_vec2=eigen_vec[match(pheno2$IID,eigen_vec[,1]),1:(no.pc+1)]

covar2=covar[,-1]
covar2=cbind(covar2,eigen_vec2[,-1])
covar2$year_sex=covar2$Year_of_birth*covar2$Sex

pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
covar2=covar2[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]
phenonames=dir("/pl/active/KellerLab/Yongkang/ldsc_results")
phenonames=phenonames[-which(phenonames%in%c("folders.csv","id_collection.csv","subtract_compare.csv"))]
phenoname=phenonames[idx]

for(chr in 1:22){
  dat=fread(paste("/rc_scratch/yoki5348/UK_PC_each_chr/not_chr",chr,".eigenvec",sep=""),header=T,stringsAsFactors=FALSE)  #Eigen vector from flash PCA 
  dat=as.data.frame(dat)
  dat=dat[match(pheno2$IID,dat$IID),]
  assign(paste("chr",chr,"_eigen",sep=""),dat)
}

pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
covar=covar[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]

no.PC_UKB=25
covar=covar[,-(5:23)]
covar2=covar[,-1]

  if(length(dir(paste("/pl/active/KellerLab/Yongkang/ldsc_results/",phenoname,"/individual_hpgs_subtracted/",sep="")))>=1100){
    for(chr in 1:22){
      idx=1
      hpgs_file=try(fread(paste("/pl/active/KellerLab/Yongkang/ldsc_results/",phenoname,"/individual_hpgs_subtracted/chr",chr,".",idx,"_ver2.csv",sep=""),header=F,stringsAsFactors=FALSE),silent=TRUE)
      hpgs_file=as.data.frame(hpgs_file)
      for(idx in 2:50){ ### Modify needed!
        hpgs_file2=try(fread(paste("/pl/active/KellerLab/Yongkang/ldsc_results/",phenoname,"/individual_hpgs_subtracted/chr",chr,".",idx,"_ver2.csv",sep=""),header=F,stringsAsFactors=FALSE))
        hpgs_file2=as.data.frame(hpgs_file2)
        if(sum(is.na(hpgs_file2[,2]))==0){
          hpgs_file[,2:3]=hpgs_file[,2:3]+hpgs_file2[,2:3]
        }
      }
      colnames(hpgs_file)=c("IID","pgs1","pgs2")
      hpgs_file=hpgs_file[match(pheno2$IID,hpgs_file$IID),]
      covar3=covar2
      eigen_vec2=get(paste("chr",chr,"_eigen",sep=""))[,3:(2+no.PC_UKB)]
      covar3=cbind(covar3,eigen_vec2)
      covar3$year_sex=covar3$Year_of_birth*covar3$Sex
      pgs=hpgs_file[,2]+hpgs_file[,3]
      lmfit=lm(pgs~.,data=covar3)
      fitted=predict(lmfit,newdata=covar3)
      hpgs_file[,2]=hpgs_file[,2]-fitted/2
      hpgs_file[,3]=hpgs_file[,3]-fitted/2
      assign(paste("hpgs_chr",chr,sep=""),hpgs_file)
      print(paste("Chromosome ",chr," is completed.",sep=""))
    }
    hpgs2=hpgs_chr1
    for(chr in 2:22){
      hpgs2[,2:3]=hpgs2[,2:3]+get(paste("hpgs_chr",chr,sep=""))[,2:3]
    }
    
    k=2*((var(hpgs2[,2],na.rm=TRUE)+var(hpgs2[,3],na.rm=TRUE))/2-(cov(hpgs2[,2],hpgs2[,3],use="complete.obs")))
    hpgs2[,2:3]=hpgs2[,2:3]/sqrt(k)
    
    gc_all=(var(hpgs2[,2],na.rm=TRUE)+var(hpgs2[,3],na.rm=TRUE)-1+cov(hpgs2[,2],hpgs2[,3],use="complete.obs"))/3
    write.csv(hpgs2,paste("/pl/active/KellerLab/Yongkang/ldsc_results/",phenoname,"/hpgs_pc",no.PC_UKB,"_subtracted.csv",sep=""),row.names=F,quote=F)
    #merge 1963-1972
    
    hpgs_fat=hpgs2[match(spo$IID_Fat,hpgs2$IID),]
    hpgs_mot=hpgs2[match(spo$IID_Mot,hpgs2$IID),]
    gt=(cov(hpgs_fat[,2],hpgs_mot[,2],use="complete.obs")+
          cov(hpgs_fat[,2],hpgs_mot[,3],use="complete.obs")+
          cov(hpgs_fat[,3],hpgs_mot[,2],use="complete.obs")+
          cov(hpgs_fat[,3],hpgs_mot[,3],use="complete.obs"))/4
    
    hpgs3=hpgs2[-which(hpgs2$IID%in%spo$IID_Fat|hpgs2$IID%in%spo$IID_Mot),]
    gc_w.o.spo=(var(hpgs3[,2],na.rm=TRUE)+var(hpgs3[,3],na.rm=TRUE)-1+cov(hpgs3[,2],hpgs3[,3],use="complete.obs"))/3
    
    hpgs3=hpgs2[which(hpgs2$IID%in%spo$IID_Fat|hpgs2$IID%in%spo$IID_Mot),]
    gc_spo=(var(hpgs3[,2],na.rm=TRUE)+var(hpgs3[,3],na.rm=TRUE)-1+cov(hpgs3[,2],hpgs3[,3],use="complete.obs"))/3
    results2=data.frame(PC_UKB=no.PC_UKB,gt=gt,gc_all=gc_all,gc_w.o.spo=gc_w.o.spo,gc_spo=gc_spo)
    write.csv(results2,paste("/pl/active/KellerLab/Yongkang/ldsc_results/",phenoname,"/each_chr_PC_merged_chr_subtracted.csv",sep=""),row.names=F,quote=F)
    
    
  }
  
