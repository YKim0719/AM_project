
args <- commandArgs(trailingOnly =TRUE)
idx=as.numeric(args[1])

library(data.table)
cor_inf=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_r2_cor_ver2.csv",header=T,stringsAsFactors=FALSE)
cor_inf=as.data.frame(cor_inf)
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)


#load("/pl/active/IBG/UKBiobank/PHENO/Keller/old_baskets/Keller-42049/ukb42049.RData")
#race=data.frame(IID=bd[,1],race=bd[,11774])
#british_IID=race$IID[which(race$race=="British")]
#write.csv(race,"/pl/active/KellerLab/Yongkang/PGS_web/UKB_race.csv",row.names=F,quote=F)
#pheno2=pheno2[which(pheno2$IID%in%british_IID),]
#write.table(pheno2,"/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",row.names=F,quote=F,col.names=T)
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_phenofiles_no_scaled_addmore.csv",header=T,stringsAsFactors=FALSE)

marriage_year=fread("/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/Assortment_Project/Marriage_Years/UKB_SpousePairs_MarriageYears.txt",header=T,stringsAsFactors=FALSE)
marriage_year=as.data.frame(marriage_year)
marriage_year=marriage_year[match(pheno2$IID,marriage_year$f.eid),]
marriage_year=marriage_year[,c(1,3,2,4)]
covar=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_covarfile3.txt",header=T,stringsAsFactors=FALSE)
covar=as.data.frame(covar)
covar=covar[match(pheno2$IID,covar$id),]
covar=covar[,-c(4:13)]
#spo=fread("/pl/active/KellerLab/jared/misc/AM_Meta_Analysis/Spouse_Detection/UKB_Spouses_All_Ethnicities.txt",header=T,stringsAsFactors=FALSE)
#spo=as.data.frame(spo)
#pairid=unique(spo$PairID)
#spo_pair=data.frame(IID_Fat=NA,IID_Mot=NA)
#for(id in pairid){
#  spo2=spo[which(spo$PairID==id),]
#  spo_pair2=data.frame(IID_Fat=spo2$f.eid[spo2$Sex==1],IID_Mot=spo2$f.eid[spo2$Sex==0])
#  spo_pair=rbind(spo_pair,spo_pair2)
#}
#spo_pair=spo_pair[-1,]
#write.csv(spo_pair,"/pl/active/KellerLab/Yongkang/UKB_spo_pair.csv",row.names=F,quote=F)

spo=fread("/pl/active/KellerLab/Yongkang/UKB_spo_pair.csv",header=T,stringsAsFactors=FALSE)
spo=as.data.frame(spo)
bin_range_year=3
bin_range_std=5
bin_year=TRUE; bin_std=TRUE

phenonames=colnames(pheno2)
phenoname=phenonames[idx]

#eigen_rare=fread("/pl/active/IBG/promero/rare-common-overlap/phenos_covar/rarePCs/pcs_thinnedPcs.txt",header=T,stringsAsFactors=FALSE) #Eigen vector from flash PCA with low frequency variables
#eigen_rare=as.data.frame(eigen_rare)

results=data.frame(income=NA,location=NA,PC_UKB=NA,gt=NA,gc_all=NA,gc_w.o.spo=NA,gc_spo=NA)
pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
covar=covar[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]

for(chr in 1:22){
  dat=fread(paste("/pl/active/KellerLab/Yongkang/UK_PC_each_chr/not_chr",chr,".eigenvec",sep=""),header=T,stringsAsFactors=FALSE)  #Eigen vector from flash PCA 
  dat=as.data.frame(dat)
  dat=dat[match(pheno2$IID,dat$IID),]
  assign(paste("chr",chr,"_eigen",sep=""),dat)
}
R_squared=data.frame(income=NA,location=NA,PC_UKB=NA,chr1_Rsq=NA,
                     chr2_Rsq=NA,chr3_Rsq=NA,chr4_Rsq=NA,chr5_Rsq=NA,chr6_Rsq=NA,chr7_Rsq=NA,
                     chr8_Rsq=NA,chr9_Rsq=NA,chr10_Rsq=NA,chr11_Rsq=NA,chr12_Rsq=NA,chr13_Rsq=NA,
                     chr14_Rsq=NA,chr15_Rsq=NA,chr16_Rsq=NA,chr17_Rsq=NA,chr18_Rsq=NA,chr19_Rsq=NA,
                     chr20_Rsq=NA,chr21_Rsq=NA,chr22_Rsq=NA)

for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.PC_UKB in seq(0,30,1)){
      covar2=covar[,-1]
      if(income==0){
        covar2=covar[,-22]
      }
      if(location==0){
        covar2=covar[,-c(4:21)]
      }
      r_square=c()
      
      for(chr in 1:22){
        idx=1
        
        hpgs_file=try(fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/individual_hpgs/chr",chr,".",idx,".csv",sep=""),header=T,stringsAsFactors=FALSE),silent=TRUE)
        if(class(hpgs_file)[1]!="try-error"){
          hpgs_file=as.data.frame(hpgs_file)
          for(idx in 2:50){ ### Modify needed!
            hpgs_file2=try(fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/individual_hpgs/chr",chr,".",idx,".csv",sep=""),header=T,stringsAsFactors=FALSE))
            if(class(hpgs_file2)[1]!="try-error"){
              hpgs_file2=as.data.frame(hpgs_file2)
              hpgs_file[,2:3]=hpgs_file[,2:3]+hpgs_file2[,2:3]
            }
          }
          
        }else{
          hpgs_file=try(fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/individual_hpgs/chr",chr,".",idx,"_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE),silent=TRUE)
          if(class(hpgs_file)[1]!="try-error"){
            hpgs_file=as.data.frame(hpgs_file)
            for(idx in 2:50){
              hpgs_file2=try(fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/individual_hpgs/chr",chr,".",idx,"_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE))
              if(class(hpgs_file2)[1]!="try-error"){
                hpgs_file2=as.data.frame(hpgs_file2)
                hpgs_file[,2:3]=hpgs_file[,2:3]+hpgs_file2[,2:3]
              }
            }
            
          }else{
            hpgs_file=try(fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/individual_hpgs/chr",chr,".",idx,"_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE),silent=TRUE)
            if(class(hpgs_file)[1]!="try-error"){
              hpgs_file=as.data.frame(hpgs_file)
              for(idx in 2:50){
                hpgs_file2=try(fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/individual_hpgs/chr",chr,".",idx,"_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE))
                if(class(hpgs_file2)[1]!="try-error"){
                hpgs_file2=as.data.frame(hpgs_file2)
                hpgs_file[,2:3]=hpgs_file[,2:3]+hpgs_file2[,2:3]
                }
              }
              
          }
          
        }
        }
        
         hpgs_file=hpgs_file[match(pheno2$IID,hpgs_file$IID),]
        covar3=covar2
        if(no.PC_UKB!=0){
          eigen_vec2=get(paste("chr",chr,"_eigen",sep=""))[,3:(2+no.PC_UKB)]
          covar3=cbind(covar3,eigen_vec2)
        }
        covar3$year_sex=covar3$Year_of_birth*covar3$Sex
        pgs=hpgs_file[,2]+hpgs_file[,3]
        lmfit=lm(pgs~.,data=covar3)
        r_square=c(r_square,summary(lmfit)$r.squared)
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
      
      #geo=fread("/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/Assortment_Project/Birthplace_Analysis/UKB_Census_Merged.txt",header=T,stringsAsFactors=FALSE)
      #geo=as.data.frame(geo)
      
      ####################
      #geo=geo[match(pheno2$IID,geo[,1]),]
      #
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
      
      k=2*((var(hpgs2[,2],na.rm=TRUE)+var(hpgs2[,3],na.rm=TRUE))/2-(cov(hpgs2[,2],hpgs2[,3],use="complete.obs")))
      hpgs2[,2:3]=hpgs2[,2:3]/sqrt(k)
      
      gc_all=(var(hpgs2[,2],na.rm=TRUE)+var(hpgs2[,3],na.rm=TRUE)-1+cov(hpgs2[,2],hpgs2[,3],use="complete.obs"))/3
      dir.create("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/")
      dir.create("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/hpgs_collection")
      
      dir.create(paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/hpgs_collection/",phenoname,sep=""))
      write.csv(hpgs2,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/hpgs_collection/",phenoname,"/hpgs_pc",no.PC_UKB,".csv",sep=""),row.names=F,quote=F)
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
      results2=data.frame(income=income,location=location,PC_UKB=no.PC_UKB,gt=gt,gc_all=gc_all,gc_w.o.spo=gc_w.o.spo,gc_spo=gc_spo)
      results=rbind(results,results2)
      R_squared2=c(income,location,no.PC_UKB,r_square)
      names(R_squared2)=c("income","location","PC_UKB",paste("chr",1:22,"_Rsq",sep=""))
      R_squared=rbind(R_squared,R_squared2)
      write.csv(results,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/gc_gt_compare_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)
      write.csv(R_squared,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/R_square_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)
    }
  }
}
