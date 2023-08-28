args <- commandArgs(trailingOnly =TRUE)
idx=as.numeric(args[1])

###Cleaning Codes to make reproducible results

##singularity run /projects/lessem/singularity/openmx.sif Rscript /pl/active/KellerLab/Yongkang/AM_source_code/find_SEM_results_include_year_sex.R $idx

#####Generate the new external summary statistic bassed results ####
library(data.table)
source("/pl/active/KellerLab/Yongkang/AM_source_code/source_code_for_SEM_ver2.R")
source("/pl/active/KellerLab/Yongkang/AM_source_code/source_code_for_various_hpgs_handle_ver2.R")
cor_inf=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_r2_cor_ver2.csv",header=T,stringsAsFactors=FALSE)
cor_inf=as.data.frame(cor_inf)
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british_include_new_alz.txt",header=T,stringsAsFactors=FALSE)
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
eigen_vec2=eigen_vec[match(pheno2$IID,eigen_vec[,1]),1:(no.pc+1)]

covar2=covar[,-1]
covar2=cbind(covar2,eigen_vec2[,-1])
covar2$year_sex=covar2$Year_of_birth*covar2$Sex

pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
covar2=covar2[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]

phenonames=colnames(pheno2)
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

for(chr in 1:22){
  dat=fread(paste("/rc_scratch/yoki5348/UK_PC_each_chr/not_chr",chr,".eigenvec",sep=""),header=T,stringsAsFactors=FALSE)  #Eigen vector from flash PCA 
  dat=as.data.frame(dat)
  dat=dat[match(pheno2$IID,dat$IID),]
  assign(paste("chr",chr,"_eigen",sep=""),dat)
}


pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
covar=covar[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]

for(idx in 6:ncol(pheno2)){
  phenoname=colnames(pheno2)[idx]
  results=data.frame(income=NA,location=NA,PC_UKB=NA,gt=NA,gc_all=NA,gc_w.o.spo=NA,gc_spo=NA)
  
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
        dir.create("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_new_spouse/")
        dir.create("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_new_spouse/hpgs_collection")
        
        dir.create(paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_new_spouse/hpgs_collection/",phenoname,sep=""))
        write.csv(hpgs2,paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_new_spouse/hpgs_collection/",phenoname,"/hpgs_pc",no.PC_UKB,".csv",sep=""),row.names=F,quote=F)
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
        write.csv(results,paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_new_spouse/gc_gt_compare_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)
        write.csv(R_squared,paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_new_spouse/R_square_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)
      }
    }
  }
  
}
#####################Plot results ################################
library(stringr)
aa=dir("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/results_external_summary/optimal_PC/gc_gt_new_spouse/")
aa=str_replace(aa,"R_square_","")
aa=str_replace(aa,"_each_chr_PC_merged_chr.csv","")
aa=str_replace(aa,"gc_gt_compare_","")
aa=aa[1:26]
idx=1
dat=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/results_external_summary/optimal_PC/gc_gt_new_spouse/gc_gt_compare_",aa[idx],"_each_chr_PC_merged_chr.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat=dat[-1,]
plot(dat$PC_UKB,dat$gt,pch=16,)
no.PC_UKB=25
for(chr in 1:22){
  idx=1
  hpgs_file=fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/regenie_train_test/pheno",phenoname,"/hpgs_results/chr",chr,".idx",idx,".csv",sep=""),header=T,stringsAsFactors=FALSE)
  hpgs_file=as.data.frame(hpgs_file)
  for(idx in 2:50){
    hpgs_file2=fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/regenie_train_test/pheno",phenoname,"/hpgs_results/chr",chr,".idx",idx,".csv",sep=""),header=T,stringsAsFactors=FALSE)
    hpgs_file2=as.data.frame(hpgs_file2)
    hpgs_file[,2:3]=hpgs_file[,2:3]+hpgs_file2[,2:3]
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
 # r_square=c(r_square,summary(lmfit)$r.squared)
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

phenonames=colnames(pheno2)
idx=13
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

k=2*((var(hpgs2[,2],na.rm=TRUE)+var(hpgs2[,3],na.rm=TRUE))/2-(cov(hpgs2[,2],hpgs2[,3],use="complete.obs")))
hpgs2[,2:3]=hpgs2[,2:3]/sqrt(k)

gc_all=(var(hpgs2[,2],na.rm=TRUE)+var(hpgs2[,3],na.rm=TRUE)-1+cov(hpgs2[,2],hpgs2[,3],use="complete.obs"))/3
#merge 1963-1972

hpgs_fat=hpgs2[match(spo$IID_fat,hpgs2$IID),]
hpgs_mot=hpgs2[match(spo$IID_mot,hpgs2$IID),]
gt=(cov(hpgs_fat[,2],hpgs_mot[,2],use="complete.obs")+
      cov(hpgs_fat[,2],hpgs_mot[,3],use="complete.obs")+
      cov(hpgs_fat[,3],hpgs_mot[,2],use="complete.obs")+
      cov(hpgs_fat[,3],hpgs_mot[,3],use="complete.obs"))/4

hpgs3=hpgs2[-which(hpgs2$IID%in%spo$IID_fat|hpgs2$IID%in%spo$IID_mot),]
gc_w.o.spo=(var(hpgs3[,2],na.rm=TRUE)+var(hpgs3[,3],na.rm=TRUE)-1+cov(hpgs3[,2],hpgs3[,3],use="complete.obs"))/3

hpgs3=hpgs2[which(hpgs2$IID%in%spo$IID_fat|hpgs2$IID%in%spo$IID_mot),]
gc_spo=(var(hpgs3[,2],na.rm=TRUE)+var(hpgs3[,3],na.rm=TRUE)-1+cov(hpgs3[,2],hpgs3[,3],use="complete.obs"))/3
results2=data.frame(income=income,location=location,PC_UKB=no.PC_UKB,gt=gt,gc_all=gc_all,gc_w.o.spo=gc_w.o.spo,gc_spo=gc_spo)
results=rbind(results,results2)
R_squared2=c(income,location,no.PC_UKB,r_square)
names(R_squared2)=c("income","location","PC_UKB",paste("chr",1:22,"_Rsq",sep=""))
R_squared=rbind(R_squared,R_squared2)
write.csv(results,paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_compare_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)
write.csv(R_squared,paste("/rc_scratch/yoki5348/UK_PC_each_chr/R_square_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)


exclude=which(is.na(hpgs2$pgs1)|is.na(pheno2[,phenoname])|is.na(pheno_dat[,2]))
if(length(exclude)!=0){
  hpgs2=hpgs2[-exclude,]
  covar2=covar2[-exclude,]
  marriage_year=marriage_year[-exclude,]
  pheno2=pheno2[-exclude,]
  pheno_dat=pheno_dat[-exclude,]
}


year_check=data.frame(years=as.numeric(names(table(marriage_year$Participant.Marriage.Year))),no.samp=as.numeric(table(marriage_year$Participant.Marriage.Year)))
year_check=cbind(year_check,year_check$years)
colnames(year_check)[3]="assigned_year"
year_check[,3]=as.character(year_check[,3])

year_check$assigned_year[year_check$years<=1972]="1963-1972"
for(id in seq(1973,1999,3)){
  year_check$assigned_year[which(year_check$years>=id&year_check$years<id+3)]=paste(id,"-",id+2,sep="")
}
year_check$plot_year=NA
aa=unique(year_check$assigned_year)
for(idx in 1:length(aa)){
  bb=year_check$years[year_check$assigned_year==aa[idx]]
  year_check$plot_year[year_check$assigned_year==aa[idx]]=round(mean(bb),2)
}
year_check_parents=data.frame(years=as.numeric(names(table(marriage_year$Parents.Marriage.Year))),no.samp=as.numeric(table(marriage_year$Parents.Marriage.Year)))
year_check_parents=cbind(year_check_parents,year_check_parents$years)
colnames(year_check_parents)[3]="assigned_year"
year_check_parents[,3]=as.character(year_check_parents[,3])
year_check_parents$assigned_year[year_check_parents$years>=1963&year_check_parents$years<=1972]="1963-1972"
year_check_parents$assigned_year[year_check_parents$years<=1935]="<1936"
for(id in seq(1936,1962,3)){
  year_check_parents$assigned_year[which(year_check_parents$years>=id&year_check_parents$years<id+3)]=paste(id,"-",id+2,sep="")
}

year_check_parents$plot_year=NA
aa=unique(year_check_parents$assigned_year)
for(idx in 1:length(aa)){
  bb=year_check_parents$years[year_check_parents$assigned_year==aa[idx]]
  year_check_parents$plot_year[year_check_parents$assigned_year==aa[idx]]=round(mean(bb),2)
}

marriage_year2=marriage_year
for(idx in 1:nrow(year_check)){
  marriage_year2$Participant.Marriage.Year[which(marriage_year$Participant.Marriage.Year==year_check$year[idx])]=year_check$plot_year[idx]
}
for(idx in 1:nrow(year_check_parents)){
  marriage_year2$Parents.Marriage.Year[which(marriage_year$Parents.Marriage.Year==year_check_parents$year[idx])]=year_check_parents$plot_year[idx]
}

aa=unique(year_check$plot_year)
year_check_gt=data.frame(years=aa[1],no.samp=sum(year_check$no.samp[year_check$plot_year==aa[1]]),assigned_year=year_check$assigned_year[which(year_check$plot_year==aa[1])[1]])
for(year in aa[-1]){
  year_check_gt2=data.frame(years=year,no.samp=sum(year_check$no.samp[year_check$plot_year==year]),assigned_year=year_check$assigned_year[which(year_check$plot_year==year)[1]])
  year_check_gt=rbind(year_check_gt,year_check_gt2)
}

aa=unique(year_check_parents$plot_year)
year_check_gc=data.frame(years=aa[1],no.samp=sum(year_check_parents$no.samp[year_check_parents$plot_year==aa[1]]),assigned_year=year_check_parents$assigned_year[which(year_check_parents$plot_year==aa[1])[1]])
for(year in aa[-1]){
  year_check_gc2=data.frame(years=year,no.samp=sum(year_check_parents$no.samp[year_check_parents$plot_year==year]),assigned_year=year_check_parents$assigned_year[which(year_check_parents$plot_year==year)[1]])
  year_check_gc=rbind(year_check_gc,year_check_gc2)
}
year_check_gt
marriage_year_fat=marriage_year2[match(spo$IID_Fat,marriage_year2$f.eid),]
marriage_year_mot=marriage_year2[match(spo$IID_Mot,marriage_year2$f.eid),]
dir_name=paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_new_spouse/SEM_analysis_results/PC25",sep="")    
dir_name2=paste("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_new_spouse/SEM_analysis_results/PC25_Rdata",sep="")    
dir.create(dir_name)
dir.create(dir_name2)
spo_fat=NA;spo_mot=NA;bin_idx=NA

pheno_fat=pheno_dat[match(spo$IID_fat,pheno_dat$IID),2]
pheno_mot=pheno_dat[match(spo$IID_mot,pheno_dat$IID),2]

if(phenoname%in%c("schizophrenia","bipolar","pgs_major_depression_non_UKB","anxiety_cc","ocd","PTSD","alzheimer")){
  g_result=try(perform_SEM_gc_gt_same_with_conti_w.o.y_only_spo(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  global_g_result=try(perform_SEM_global_with_conti_w.o.y_only_spo(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  trend_result=try(perform_SEM_linear_trend_with_conti_w.o.y_only_spo(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  nonlinear_trend_result=try(perform_SEM_nonlinear_trend_with_conti_w.o.y_only_spo(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  cor_ref=cor_inf$cor_ref[cor_inf$pheno==phenoname]
  Omega2mu_gc_gt=mxEval(4*gt1/se_ref-cor_ref,global_g_result)
  Omega2mu_gc_gt.se=mxSE(4*gt1/se_ref-cor_ref,global_g_result)
  gtest_1967=mxEval(2*gc11-gc1,nonlinear_trend_result)
  gtest_1967.se=mxSE(2*gc11-gc1,nonlinear_trend_result)
  gcest_1937=mxEval((gt1+gc11)/2,nonlinear_trend_result)
  gcest_1937.se= mxSE((gt1+gc11)/2,nonlinear_trend_result)
  rhoest_1967=mxEval(2*gcest11-gcest1,nonlinear_trend_result)
  rhoest_1967.se=mxSE(2*gcest11-gcest1,nonlinear_trend_result)
  rhoest_1937=mxEval((gtest1+gcest11)/2,nonlinear_trend_result)
  rhoest_1937.se= mxSE((gtest1+gcest11)/2,nonlinear_trend_result)
}else{
  g_result=try(perform_SEM_gc_gt_same_with_conti_only_spo(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  global_g_result=try(perform_SEM_global_with_conti_only_spo(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  trend_result=try(perform_SEM_linear_trend_with_conti_only_spo(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  nonlinear_trend_result=try(perform_SEM_nonlinear_trend_with_conti_only_spo(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  Omega2mu_gc_gt=mxEval(gt1-Omega1^2*spouse,global_g_result)
  Omega2mu_gc_gt.se=mxSE(gt1-Omega1^2*spouse,global_g_result)
  gtest_1967=mxEval(2*gc11-gc1,nonlinear_trend_result)
  gtest_1967.se=mxSE(2*gc11-gc1,nonlinear_trend_result)
  gcest_1937=mxEval((gt1+gc11)/2,nonlinear_trend_result)
  gcest_1937.se= mxSE((gt1+gc11)/2,nonlinear_trend_result)
  R2_UKB=cor_inf$R2_UKB[cor_inf$pheno==phenoname]
  rhoest_1967=mxEval(4*(2*gc11-gc1)/R2_UKB,nonlinear_trend_result)
  rhoest_1967.se=mxSE(4*(2*gc11-gc1)/R2_UKB,nonlinear_trend_result)
  rhoest_1937=mxEval(4*(gt1+gc11)/(2*R2_UKB),nonlinear_trend_result)
  rhoest_1937.se= mxSE(4*(gt1+gc11)/(2*R2_UKB),nonlinear_trend_result)
}

bb=mxCompare(global_g_result ,g_result)
bb2=mxCompare(nonlinear_trend_result,trend_result )
bb=rbind(bb,bb2)
write.csv(bb,paste(dir_name,"/",phenoname,"/lrt_test_results_only_spo.csv",sep=""),row.names=F,quote=F)

est_trend=data.frame(phenoname,gtest_1967,
                     gtest_1967.se,
                     gcest_1937,
                     gcest_1937.se,
                     rhoest_1967,
                     rhoest_1967,
                     rhoest_1937,
                     rhoest_1937,
                     Omega2mu_gc_gt,
                     Omega2mu_gc_gt.se
)
write.csv(est_trend,paste(dir_name,"/",phenoname,"/est_trend_only_spo.csv",sep=""),row.names=F,quote=F)
dir.create(paste(dir_name2,"/",phenoname,sep=""))
save.image(file=paste(dir_name2,"/",phenoname,"/Rimage_only_spo.RData",sep=""))

dat1=read.csv(paste(dir_name,"/",phenoname,"/diseq_AM_assumption_gc_gt_nonlinear_only_spo.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat2=read.csv(paste(dir_name,"/",phenoname,"/diseq_AM_assumption_gc_gt_linear_increasing_only_spo.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat3=read.csv(paste(dir_name,"/",phenoname,"/diseq_AM_assumption_gc_gt_global_only_spo.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat4=read.csv(paste(dir_name,"/",phenoname,"/diseq_AM_assumption_gc_gt_same_only_spo.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat5=read.csv(paste(dir_name,"/",phenoname,"/est_trend.csv",sep=""),header=T,stringsAsFactors=FALSE)

if(phenoname%in%c("schizophrenia","bipolar","pgs_major_depression_non_UKB","anxiety_cc","ocd","PTSD","alzheimer")){
  rho_gt=paste("gtest_ref",1:10,sep="")
  rho_gt_dat=dat1[match(rho_gt,dat1$name),-c(1:4)]
  rho_gc=paste("gcest_ref",1:11,sep="")
  rho_gc_dat=dat1[match(rho_gc,dat1$name),-c(1:4)]
  
}else{
  rho_gt=paste("gt",1:10,sep="")
  rho_gt_dat=dat1[match(rho_gt,dat1$matrix),-c(1:4)]*4/cor_inf$R2_UKB[cor_inf$pheno==phenoname]
  rho_gc=paste("gc",1:11,sep="")
  rho_gc_dat=dat1[match(rho_gc,dat1$matrix),-c(1:4)]*4/cor_inf$R2_UKB[cor_inf$pheno==phenoname]
  rho_gc_spo=paste("gcprime",1:11,sep="")
  rho_gc_spo_dat=dat1[match(rho_gc_spo,dat1$matrix),-c(1:4)]*4/cor_inf$R2_UKB[cor_inf$pheno==phenoname]
  
}

year_gt=c("1963-1972",paste(seq(1973,1999,3),seq(1975,2001,3),sep="-"))
year_gc=c("<1936",paste(seq(1936,1962,3),seq(1938,1964,3),sep="-"),"1963-1972")

year_line1=c("<1936",paste(seq(1936,1962,3),"-",sep=""),"1963-",paste(seq(1973,1999,3),"-",sep=""))
year_line2=c("",seq(1938,1964,3),"1972",seq(1975,2001,3))
x_points=1:20

####Draw g plot along the marriage year ###########
dir.create(paste(dir_name,"/",phenoname,sep=""))


png(paste(dir_name,"/",phenoname,"/rhoplot_fiex_axis_only_spo.png",sep=""),width=1500,height=800,pointsize=16)
yax=c(rho_gt_dat$lower.c.i,rho_gt_dat$upper.c.i,rho_gc_dat$lower.c.i,rho_gc_dat$upper.c.i,0,cor_inf$cor_ref[which(cor_inf$pheno==phenoname)],cor_inf$cor_UKB[which(cor_inf$pheno==phenoname)])
yax=na.omit(yax)
plot(11:20,rho_gt_dat$Estimate,pch=16,type="p",xlab="Marriage year",col="red",ylab="rho",main=phenoname,xlim=c(1,20),ylim=c(min(yax),max(yax)),axes=FALSE,cex=2)
lines(11:20,rho_gt_dat$Estimate,col="red",lwd=2)
lines(11:20,rho_gt_dat$lower.c.i,lty=2,col="red",lwd=2)
lines(11:20,rho_gt_dat$upper.c.i,lty=2,col="red",lwd=2)
if(phenoname%in%c("schizophrenia","bipolar","pgs_major_depression_non_UKB","anxiety_cc","ocd","PTSD","alzheimer")){
  lines(11:20,c(dat2[dat2$matrix=="gt1",5]+(0:9)*dat2[dat2$matrix=="bgt",5])*4/cor_inf$R2_ref[cor_inf$pheno==phenoname],col="brown",lwd=2.5)
}else{
  lines(11:20,c(dat2[dat2$matrix=="gt1",5]+(0:9)*dat2[dat2$matrix=="bgt",5])*4/cor_inf$R2_UKB[cor_inf$pheno==phenoname],col="brown",lwd=2.5)
  
}

abline(a=0,b=0,col="black",lty=2,lwd=1.5)
axis(1,1:20,year_line1,line=0)
axis(1,1:20,year_line2,line=1,tick=FALSE)
axis(2)
box()

points(1:11,rho_gc_dat$Estimate,pch=16,col="skyblue",cex=2)
lines(1:11,rho_gc_dat$Estimate,pch=16,col="skyblue",lwd=2)
lines(1:11,rho_gc_dat$upper.c.i,col="skyblue",lwd=2,lty=2)
lines(1:11,rho_gc_dat$lower.c.i,col="skyblue",lwd=2,lty=2)
if(phenoname%in%c("schizophrenia","bipolar","pgs_major_depression_non_UKB","anxiety_cc","ocd","PTSD","alzheimer")){
  lines(1:11,c(dat2[dat2$matrix=="gc1",5]+(0:10)*dat2[dat2$matrix=="bgc",5])*4/cor_inf$R2_ref[cor_inf$pheno==phenoname],col="#008B8B",lwd=2.5)
}else{
  lines(1:11,c(dat2[dat2$matrix=="gc1",5]+(0:10)*dat2[dat2$matrix=="bgc",5])*4/cor_inf$R2_UKB[cor_inf$pheno==phenoname],col="#008B8B",lwd=2.5)
  
}

if(phenoname%in%c("schizophrenia","bipolar","pgs_major_depression_non_UKB","anxiety_cc","ocd","PTSD","alzheimer")){
  abline(a=cor_inf$cor_ref[which(cor_inf$pheno==phenoname)],b=0,col="green",lwd=2.5,lty=2)
}else{
  abline(a=cor_inf$cor_UKB[which(cor_inf$pheno==phenoname)],b=0,col="green",lwd=2.5,lty=2)
  
}

dev.off()




#####Generate the new results for comparing internal analysis and that of external
#Generate rg values
library(data.table)
library(stringr)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)


check_count=data.frame(phenoname=NA,count=NA)

for(idx in 7:ncol(pheno2)){
  phenoname=colnames(pheno2)[idx]
  check_count2=data.frame(phenoname=phenoname,count=  sum(str_detect(dir(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,sep="")),paste(".regenie",sep="")))
  )
  
  check_count=rbind(check_count,check_count2)
}

for(idx in which(colnames(pheno2)%in%check_count[check_count$count==22,1])){
  phenoname=colnames(pheno2)[idx]
  chr=1
  dat=fread(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoid,"/test_chr",chr,"_",phenoname,".regenie",sep=""),header=T,stringsAsFactors=FALSE)
  dat=as.data.frame(dat)
  for(chr in 2:22){
    dat2=fread(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoid,"/test_chr",chr,"_",phenoname,".regenie",sep=""),header=T,stringsAsFactors=FALSE)
    dat2=as.data.frame(dat2)
    dat=rbind(dat,dat2)
  }
  results=data.frame(snpid=dat[,3],hg18chr=dat[,1],bp=dat[,2],
                     a1=dat[,5],a2=dat[,4],beta=dat$BETA,se=dat$SE,pval=10^(-dat$LOG10P))
  n=mean(dat$N)
  write.table(results,paste( "/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/sumstat.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  
  #check_where=data.frame(phenoname=NA,general_loc=NA)
  #for(idx in 7:ncol(pheno2)){
  #  phenoname=colnames(pheno2)[idx]
  #  check_where2=data.frame(phenoname=phenoname,general_loc=sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))=="sumstat_pheno_all.txt")==1)
  #  check_where=rbind(check_where,check_where2)
  #}
  system(paste("singularity run /pl/active/KellerLab/Yongkang/ldsc3.sif  munge_sumstats.py ",
               "--sumstats /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/sumstat.txt --N ",n,
               " --out /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc",sep=""))
  if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))=="sumstat_pheno_all.txt")==1){
    dat=fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/sumstat_pheno_all.txt",sep=""),header=T,stringsAsFactors=FALSE)
    dat=as.data.frame(dat)
    system(paste("singularity run /pl/active/KellerLab/Yongkang/ldsc3.sif  munge_sumstats.py ",
                 "--sumstats /pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/sumstat_pheno_all.txt --N ",200000,
                 " --out /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc_external",sep=""))
    system(paste("singularity run /pl/active/KellerLab/Yongkang/ldsc3.sif ldsc.py --rg /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc.sumstats.gz,/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc_external.sumstats.gz --ref-ld-chr /pl/active/KellerLab/Yongkang/hapmap3/ldcomp/chr@ --w-ld-chr /pl/active/KellerLab/Yongkang/hapmap3/ldcomp/chr@ --out /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc",sep=""))
    
    
  }
  if(sum(phenoname %in% c("AgeOfInitiation_gscan","CigarettesPerDay_gscan","SmokingInitiation_gscan","SmokingCessation_gscan","SmokingCessation_gscan","DrinksPerWeek_gscan"))==1){
    phenoname2=str_replace(phenoname,"_gscan","")
    dat=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/prs_external/summary_statistics/gscan/",phenoname2,".WithoutUKB.txt",sep=""),header=T,stringsAsFactors=FALSE)
    dat=as.data.frame(dat)
    dat_merged_ref2=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/matched_snp_info.txt",header=T,stringsAsFactors=FALSE)
    dat_merged_ref2=as.data.frame(dat_merged_ref2)
    dat_matched=dat[match(dat_merged_ref2$rsid,dat$RSID),]
    dat_matched=cbind(dat_matched,dat_merged_ref2[,c(3,5)])
    info_snp2=dat_matched[,c("SNPID_GRCH38","CHROM","REF","ALT","GRCH38_pos","N","SE","BETA","PVALUE")]
    colnames(info_snp2)=c("snpid","hg18chr","a2","a1","pos","n_eff","se","beta","pval")
    write.table(info_snp2,paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/sumstat_pheno_all.txt",sep=""),row.names=F,col.names=T,quote=F)
    n_eff=mean(info_snp2$n_eff,na.rm=TRUE)
    system(paste("singularity run /pl/active/KellerLab/Yongkang/ldsc3.sif  munge_sumstats.py ",
                 "--sumstats /pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/sumstat_pheno_all.txt --N ",n_eff,
                 " --out /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc_external",sep=""))
    system(paste("singularity run /pl/active/KellerLab/Yongkang/ldsc3.sif ldsc.py --rg /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc.sumstats.gz,/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc_external.sumstats.gz --ref-ld-chr /pl/active/KellerLab/Yongkang/hapmap3/ldcomp/chr@ --w-ld-chr /pl/active/KellerLab/Yongkang/hapmap3/ldcomp/chr@ --out /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc",sep=""))
    
  }
  
}

results=data.frame(pheno=NA,rg=NA)
for(idx in which(colnames(pheno2)%in%check_count[check_count$count==22,1])){
  phenoname=colnames(pheno2)[idx]
  aa=try(read.table(paste("/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc.log",sep=""),header=F,nrows=50,sep="!"))
  if(class(aa)[1]!="try-error"){
    results2=data.frame(pheno=phenoname,rg=aa[50,1])
    results=rbind(results,results2)
    
  }
}  
write.csv(results,"/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/rg_results_report.csv",row.names=F,quote=F)


#sinteractive --partition=atesting --qos=testing --ntasks=4  

library(data.table)
library(stringr)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)


check_count=data.frame(phenoname=NA,count=NA)
for(idx in 7:ncol(pheno2)){
  phenoname=colnames(pheno2)[idx]
  check_count2=data.frame(phenoname=phenoname,count=  sum(str_detect(dir(paste( "/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,sep="")),paste(".regenie",sep="")))
  )
  
  check_count=rbind(check_count,check_count2)
}
results=data.frame(pheno=NA,rg=NA)
for(idx in which(colnames(pheno2)%in%check_count[check_count$count>=22,1])){
  phenoname=colnames(pheno2)[idx]
  aa=try(read.table(paste("/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/results_ldsc.log",sep=""),header=F,nrows=50,sep="!"))
  if(class(aa)[1]!="try-error"){
    results2=data.frame(pheno=phenoname,rg=aa[47,1])
    results=rbind(results,results2)
    
  }
}  
write.csv(results,"/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/rg_results_report_intercept.csv",row.names=F,quote=F)



########Generate simulation results with PCs without target chromosomes

###########Remove HMC region to generate eigen vectors without molecular structural effects on chromosome 6
#sinteractive --partition=atesting --ntasks=16 --mem=30Gb --nodes=1 --qos=testing --time=03:00:00
#MHC region by hg38 chr6:28,510,120-33,480,577
#MHC region by hg37 chr6:28,477,797-33,448,354
library(data.table)
dat=fread("/pl/active/IBG/UKBiobank/GENO/QCed/genotyped/white/eur.qc.bim",header=F,stringsAsFactors=FALSE)
dat=as.data.frame(dat)
dat2=dat[dat$V1==6,]
dat2=dat2[dat2$V4>=26477797 &dat2$V4<=35448354, ]
write.table(dat2[,2],"/pl/active/KellerLab/Yongkang/HMC_region_variant_orig_UKB.txt",row.names=F,col.names=F,quote=F)
dim(dat2)
system("~/plink2 --bfile /pl/active/IBG/UKBiobank/GENO/QCed/genotyped/white/eur.qc --exclude /pl/active/KellerLab/Yongkang/HMC_region_variant_orig_UKB.txt --make-bed --out /pl/active/KellerLab/Yongkang/orig_UKB_exclude_HMC_region")
#######Perform deconvolution to remove gc inflation
~/plink2 --bfile /pl/active/IBG/UKBiobank/GENO/QCed/genotyped/white/eur.qc --threads 16 --exclude /pl/active/KellerLab/Yongkang/HMC_region_variant_orig_UKB.txt --make-bed --out /pl/active/KellerLab/Yongkang/orig_UKB_exclude_HMC_region


##########Check relationship between PCs and geographical information & income
##Need to check PCs
library(data.table)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)

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
dir.create("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/SEM_analysis_results/")
eigen_vec2=eigen_vec[match(pheno2$IID,eigen_vec[,1]),1:(no.pc+1)]

covar2=covar[,-1]
covar2=cbind(covar2,eigen_vec2[,-1])
covar2$year_sex=covar2$Year_of_birth*covar2$Sex

##########Perform Internal analysis#######################

library(data.table)
library(stringr)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)

###########Check not running #########################
check_count_ver2=data.frame(phenoname=NA,chr=NA,count=NA,phenoidx=NA)

for(idx in 7:ncol(pheno2)){
  phenoname=colnames(pheno2)[idx]
  for(chr in 1:22){
    check_count2=data.frame(phenoname=phenoname,chr=chr,count=  sum(str_detect(dir(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,sep="")),paste("test_chr",chr,"_",phenoname,".regenie",sep=""))),
                            phenoidx=idx)
    
    check_count_ver2=rbind(check_count_ver2,check_count2)
    
  }
  
}

check_count_ver3=(check_count_ver2[check_count_ver2[,3]==0,])
check_count_ver3=check_count_ver3[-1,]


for phenoidx in {29,30,31}
do
for chr in {1,2}
do
sbatch /pl/active/KellerLab/Yongkang/AM_source_code/perform_regenie_internal_analysis_test.sh $phenoidx $chr
done
done

phenoidx=38
for chr in {1,8,9,10,11,12,13}
do
sbatch /pl/active/KellerLab/Yongkang/AM_source_code/perform_regenie_internal_analysis_test.sh $phenoidx $chr
done

#########train not run ###############
check_count_ver2=data.frame(phenoname=NA,chr=NA,count=NA,phenoidx=NA)

for(idx in 7:ncol(pheno2)){
  phenoname=colnames(pheno2)[idx]
  for(chr in 1:22){
    check_count2=data.frame(phenoname=phenoname,chr=chr,count=  sum(str_detect(dir(paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,sep="")),paste("train_chr",chr,"_",phenoname,".regenie",sep=""))),
                            phenoidx=idx)
    
    check_count_ver2=rbind(check_count_ver2,check_count2)
    
  }
  
}

check_count_ver3=(check_count_ver2[check_count_ver2[,3]==0,])
check_count_ver3=check_count_ver3[-1,]
chr=2
for phenoidx in {20,27}
do
sbatch /pl/active/KellerLab/Yongkang/AM_source_code/perform_regenie_internal_analysis_train.sh $phenoidx $chr
done
##############################################

dir.create("/pl/active/KellerLab/Yongkang/hapmap3/ldfile")
NCORES <- 20
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


library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)
library(stringr)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)
phenoname=colnames(pheno2)[idx]

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
mapfile=fread("/pl/active/KellerLab/Yongkang/hapmap3/bed_files/ldfile_info.txt",header=T,stringsAsFactors=FALSE)
mapfile=as.data.frame(mapfile)
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

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
dat_merged$beta_ldpred_inf=beta_inf
fwrite(dat_merged,paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"ldpred_results.csv",sep=""),row.names=F,quote=F)

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
fwrite(dat_merged,paste( "/pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"ldpred_results.csv",sep=""),row.names=F,quote=F)

system(paste("rm ",tmp, ".sbk",sep=""))
system(paste("rm /pl/active/KellerLab/Yongkang/PGS_web/regenie_internal_pheno_all/",phenoname,"/corr_file_final_bipolar_file1.RData",sep=""))
###Check height results from hapmap3 only results #########
library(data.table)
library(stringr)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)

chr=1
dat_hapmap3=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/hapmap3_chr",chr,".txt",sep=""),header=T,stringsAsFactors=FALSE)
dat_hapmap3=as.data.frame(dat_hapmap3)
for(chr in 2:22){
  dat2=fread(paste("/pl/active/KellerLab/Yongkang/UKBhap/Meng/GRCH38_37_matched_info/hapmap3_chr",chr,".txt",sep=""),header=T,stringsAsFactors=FALSE)
  dat2=as.data.frame(dat2)
  dat_hapmap3=rbind(dat_hapmap3,dat2)

}


phenoname="height"
chr=1
dat_internal=fread(paste( "/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/chr",chr,"_",phenoname,".regenie",sep=""),header=T,stringsAsFactors=FALSE)
dat_internal=as.data.frame(dat_internal)
for(chr in 2:22){
  dat2=fread(paste( "/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/",phenoname,"/chr",chr,"_",phenoname,".regenie",sep=""),header=T,stringsAsFactors=FALSE)
  dat2=as.data.frame(dat2)
  dat_internal=rbind(dat_internal,dat2)
}
dim(dat_internal)#1308282
dat_internal=dat_internal[match(dat_hapmap3$SNPID_GRCH38,dat_internal$ID),]
sum(!is.na(dat_internal$ID))#1038396
results=data.frame(snpid=dat_internal[,3],hg18chr=dat_internal[,1],bp=dat_internal[,2],
                   a1=dat_internal[,5],a2=dat_internal[,4],beta=dat_internal$BETA,se=dat_internal$SE,pval=10^(-dat_internal$LOG10P))
results=na.omit(results)
n=mean(dat_internal$N,na.rm=TRUE)
dir.create("/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3")
write.table(results,paste( "/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/sumstat.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

system(paste("singularity run /pl/active/KellerLab/Yongkang/ldsc3.sif  munge_sumstats.py ",
             "--sumstats /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/sumstat.txt --N ",n,
             " --out /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/results_ldsc",sep=""))

if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))=="sumstat_pheno_all.txt")==1){
  dat_external=fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/sumstat_pheno_all.txt",sep=""),header=T,stringsAsFactors=FALSE)
  dat_external=as.data.frame(dat_external)
  dim(dat_external)#1772992
  dat_external=dat_external[match(dat_hapmap3$SNPID_GRCH38,dat_external$snpid),]
  
  sum(!is.na(dat_external$snpid))#923007
  dat_external=na.omit(dat_external)
  write.table(dat_external,"/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/sumstat_external.txt",row.names=F,col.names=T,quote=F)
  system(paste("singularity run /pl/active/KellerLab/Yongkang/ldsc3.sif  munge_sumstats.py ",
               "--sumstats /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/sumstat_external.txt --N ",200000,
               " --out /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/results_ldsc_external",sep=""))
  system(paste("singularity run /pl/active/KellerLab/Yongkang/ldsc3.sif ldsc.py --rg /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/results_ldsc.sumstats.gz,/pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/results_ldsc_external.sumstats.gz --ref-ld-chr /pl/active/KellerLab/Yongkang/hapmap3/ldcomp/chr@ --w-ld-chr /pl/active/KellerLab/Yongkang/hapmap3/ldcomp/chr@ --out /pl/active/KellerLab/Yongkang/results_LDSC_compare_sumstats/height_only_hapmap3/results_ldsc_rg",sep=""))
  
  
}


##Update slide with R^2 values
##Perform ld-pred2 for discovery samples
############Generate hPGS based on the internal analysis


