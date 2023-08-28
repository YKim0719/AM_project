library(data.table)
library(stringr)
aa=dir("/pl/active/IBG/hrc_ega/EGAD00001002729")[str_detect(dir("/pl/active/IBG/hrc_ega/EGAD00001002729"),".legend" )]
dat2=fread("/pl/active/IBG/UKBiobank/GENO/QCed/genotyped/white/eur.qc.bim",header=F,stringsAsFactors=FALSE)
dat2=as.data.frame(dat2)
write.table(dat2$V2,"/pl/active/KellerLab/Yongkang/HRC_pca/snp_list_UKB.txt",row.names=F,col.names=F,quote=F)
for(chr in 2:22){
  dat=fread(paste("/pl/active/IBG/hrc_ega/EGAD00001002729/",aa[chr],sep=""),header=T,stringsAsFactors=FALSE)
  dat=as.data.frame(dat)
  dat3=dat2[dat2$V1==chr,]
  
  dat4=dat[dat$position%in%dat3$V4,]
  
  dat5=dat3[match(dat4$position,dat3$V4),]
  
  snp_check=c()
  for(idx in 1:nrow(dat4)){
    snp_check=c(snp_check,as.numeric(sum(c(dat4$a0[idx],dat4$a1[idx])%in%c(dat5$V5[idx],dat5$V6[idx]))==2))
  }
  dat_used=dat4[which(snp_check==1),]
  dir.create("/pl/active/KellerLab/Yongkang/HRC_pca")
  write.table(dat_used$id,paste("/pl/active/KellerLab/Yongkang/HRC_pca/snp_list_chr",chr,".txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  write.table(dat_used,paste("/pl/active/KellerLab/Yongkang/HRC_pca/snp_info_chr",chr,".txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  
}

aa=dir("/pl/active/IBG/hrc_ega/EGAD00001002729")[str_detect(dir("/pl/active/IBG/hrc_ega/EGAD00001002729"),".vcf" )]
aa=aa[-which(str_detect(aa,"tbi"))]

chr=1
for(chr in 2:22){
  system(paste("~/plink2 --vcf /pl/active/IBG/hrc_ega/EGAD00001002729/",aa[chr]," --make-bed --extract ","/pl/active/KellerLab/Yongkang/HRC_pca/snp_list_UKB.txt --out ","/pl/active/KellerLab/Yongkang/HRC_pca/chr",chr,sep=""))
  
}

chr=1
dat=fread(paste("/pl/active/KellerLab/Yongkang/HRC_pca/chr",chr,".bim",sep=""),header=F,stringsAsFactors=FALSE)
dat=as.data.frame(dat)
exclude=as.numeric(names(table(dat$V4))[which(table(dat$V4)>=2)])
lists=unique(dat$V2[dat$V4%in%exclude])
for(chr in 2:22){
  dat=fread(paste("/pl/active/KellerLab/Yongkang/HRC_pca/chr",chr,".bim",sep=""),header=F,stringsAsFactors=FALSE)
  dat=as.data.frame(dat)
  exclude=as.numeric(names(table(dat$V4))[which(table(dat$V4)>=2)])
  lists=c(lists,unique(dat$V2[dat$V4%in%exclude]))
  
}
write.table(lists,"/pl/active/KellerLab/Yongkang/HRC_pca/excluded.txt",row.names=F,quote=F)
aa=data.frame(paste("/pl/active/KellerLab/Yongkang/HRC_pca/chr",2:22,".bed",sep=""),
           paste("/pl/active/KellerLab/Yongkang/HRC_pca/chr",2:22,".bim",sep=""),
           paste("/pl/active/KellerLab/Yongkang/HRC_pca/chr",2:22,".fam",sep=""))
write.table(aa,"/pl/active/KellerLab/Yongkang/HRC_pca/merge_list.txt",row.names=F,col.names=F,quote=F)


for(chr in 1:22){
  system(paste("~/plink2 --bfile /pl/active/KellerLab/Yongkang/HRC_pca/chr",chr," --exclude /pl/active/KellerLab/Yongkang/HRC_pca/excluded.txt --make-bed --out /pl/active/KellerLab/Yongkang/HRC_pca/chr",chr,sep=""))
  
}
system(paste("~/plink --noweb --bfile /pl/active/KellerLab/Yongkang/HRC_pca/chr1 --merge-list /pl/active/KellerLab/Yongkang/HRC_pca/merge_list.txt --make-bed --out /pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca"))
##########Remove variants which have Ambiguous strands##############
system(paste("~/plink2 --bfile /pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca --freq  --out /pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca"))
system(paste("~/plink2 --bfile /pl/active/IBG/UKBiobank/GENO/QCed/genotyped/white/eur.qc --freq  --out /pl/active/KellerLab/Yongkang/HRC_pca/UKB"))
dat=fread("/pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca.afreq",header=T,stringsAsFactors=FALSE)
dat=as.data.frame(dat)
dat2=fread("/pl/active/KellerLab/Yongkang/HRC_pca/UKB.afreq",header=T,stringsAsFactors=FALSE)
dat2=as.data.frame(dat2)
dat2=dat2[match(dat$ID,dat2$ID),]

if(sum(dat$ALT!=dat2$ALT)!=0){
  dat2$ALT_FREQS[which(dat$ALT!=dat2$ALT)]=1-dat2$ALT_FREQS[which(dat$ALT!=dat2$ALT)]
}
if(sum(abs(dat$ALT_FREQS-dat2$ALT_FREQS)>=0.1)!=0){
  remove_list=dat2$ID[which(abs(dat$ALT_FREQS-dat2$ALT_FREQS)>=0.1)]
}
write.table(remove_list,"/pl/active/KellerLab/Yongkang/HRC_pca/remove_list_ambiuous_strands.txt",row.names=F,quote=F,col.names=F)
############################
system(paste("~/plink2 --nonfounders --exclude /pl/active/KellerLab/Yongkang/HRC_pca/remove_list_ambiuous_strands.txt --bfile /pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca --maf 0.05 --freq counts --pca 50 allele-wts --out /pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca",sep=""))

system(paste("~/plink2 --nonfounders --bfile /pl/active/IBG/UKBiobank/GENO/QCed/genotyped/white/eur.qc",
             " --read-freq /pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca",
             ".acount --score /pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca",
             ".eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize ",
             " --score-col-nums 6-55 --out /pl/active/KellerLab/Yongkang/HRC_pca/eigen_vec_from_hrc_to_UKB_top50",sep=""))

#system(paste("~/plink2 --bfile /pl/active/KellerLab/Yongkang/HRC_pca/chr1 --pmerge-list /pl/active/KellerLab/Yongkang/HRC_pca/merge_list.txt --make-bed --out /pl/active/KellerLab/Yongkang/HRC_pca/hrc_for_pca"))
library(data.table)
cor_inf=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_r2_cor_ver2.csv",header=T,stringsAsFactors=FALSE)
cor_inf=as.data.frame(cor_inf)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)


load("/pl/active/IBG/UKBiobank/PHENO/Keller/old_baskets/Keller-42049/ukb42049.RData")
race=data.frame(IID=bd[,1],race=bd[,11774])
british_IID=race$IID[which(race$race=="British")]
pheno2=pheno2[which(pheno2$IID%in%british_IID),]
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_phenofiles_no_scaled_addmore.csv",header=T,stringsAsFactors=FALSE)

marriage_year=fread("/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/Assortment_Project/Marriage_Years/UKB_SpousePairs_MarriageYears.txt",header=T,stringsAsFactors=FALSE)
marriage_year=as.data.frame(marriage_year)
marriage_year=marriage_year[match(pheno2$IID,marriage_year$f.eid),]
marriage_year=marriage_year[,c(1,3,2,4)]
covar=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_covarfile3.txt",header=T,stringsAsFactors=FALSE)
covar=as.data.frame(covar)
covar=covar[match(pheno2$IID,covar$id),]
covar=covar[,-c(4:13)]
spo=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/PGS/210517_NEW/spouse_dat_pheno_pgs_results.csv",header=T,stringsAsFactors=FALSE)
spo=as.data.frame(spo)

bin_range_year=3
bin_range_std=5
bin_year=TRUE; bin_std=TRUE

phenoname="EA"
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/",sep=""))%in%"chr_merged_results.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/chr_merged_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))%in%"chr_merged_results.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/chr_merged_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/",sep=""))%in%"chr_merged_results_ver2.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/chr_merged_results_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))%in%"chr_merged_results_ver2.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/chr_merged_results_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
#eigen_rare=fread("/pl/active/IBG/promero/rare-common-overlap/phenos_covar/rarePCs/pcs_thinnedPcs.txt",header=T,stringsAsFactors=FALSE) #Eigen vector from flash PCA with low frequency variables
#eigen_rare=as.data.frame(eigen_rare)

eigen_hrc=fread("/pl/active/KellerLab/Yongkang/HRC_pca/eigen_vec_from_hrc_to_UKB_top50.sscore",header=T,stringsAsFactors=FALSE)
eigen_hrc=as.data.frame(eigen_hrc)

eigen_vec=fread("/pl/active/IBG/UKBiobank/exclusions/relatedness/white/pca/projections.txt",header=T,stringsAsFactors=FALSE)  #Eigen vector from flash PCA 
eigen_vec=as.data.frame(eigen_vec)

results=data.frame(income=NA,location=NA,PC_UKB=NA,PC_HRC=NA,gt=NA,gc_all=NA,gc_w.o.spo=NA,gc_spo=NA)
pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
covar=covar[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]

for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.PC_HRC in seq(0,50,5)){
      for(no.PC_UKB in seq(0,100,5)){
        if(no.PC_HRC!=0&no.PC_UKB!=0){
          eigen_hrc2=eigen_hrc[match(pheno2$IID,eigen_hrc$IID),5:(4+no.PC_HRC)]
          eigen_vec2=eigen_vec[match(pheno2$IID,eigen_vec$IID),3:(2+no.PC_UKB)]
          eigen_vec2=cbind(eigen_hrc2,eigen_vec2)  
        }
        if(no.PC_HRC!=0&no.PC_UKB==0){
          eigen_hrc2=eigen_hrc[match(pheno2$IID,eigen_hrc$IID),5:(4+no.PC_HRC)]
          eigen_vec2=eigen_hrc2
        }
        if(no.PC_HRC==0&no.PC_UKB!=0){
          eigen_vec2=eigen_vec[match(pheno2$IID,eigen_vec$IID),3:(2+no.PC_UKB)]
        }
        
        covar2=covar[,-1]
        if(income==0){
          covar2=covar[,-22]
        }
        if(location==0){
          covar2=covar[,-c(4:21)]
        }
        if(no.PC_HRC!=0|no.PC_UKB!=0){
          covar2=cbind(covar2,eigen_vec2)
        }
        covar2$year_sex=covar2$Year_of_birth*covar2$Sex
        
        
        #geo=fread("/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/Assortment_Project/Birthplace_Analysis/UKB_Census_Merged.txt",header=T,stringsAsFactors=FALSE)
        #geo=as.data.frame(geo)
        
        ####################
        #geo=geo[match(pheno2$IID,geo[,1]),]
        #
        phenonames=colnames(pheno2)
        idx=14
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
        
        hpgs2=hpgs[match(pheno2$IID,hpgs$IID),]
        pgs=c(hpgs2[,2]+hpgs2[,3])
        lmfit=lm(pgs~.,data=covar2)
        fitted=predict(lmfit,newdata=covar2)
        hpgs2[,2]=hpgs2[,2]-fitted/2
        hpgs2[,3]=hpgs2[,3]-fitted/2
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
        results2=data.frame(income=income,location=location,PC_UKB=no.PC_UKB,PC_HRC=no.PC_HRC,gt=gt,gc_all=gc_all,gc_w.o.spo=gc_w.o.spo,gc_spo=gc_spo)
        results=rbind(results,results2)
      }
    }
  }
}
results=results[-1,]
write.csv(results,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_compare_",phenoname,".csv",sep=""),row.names=F,quote=F)

###################Plot results #################################
results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_compare.csv",header=T,stringsAsFactors=FALSE)

dir.create("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB")

for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.HRC in seq(0,50,5)){
      results2=results[which(results$income==income&results$location==location&results$PC_HRC==no.HRC),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB/income",income,"_location",location,"_UKB_varying_no.HRC",no.HRC,".png",sep=""),width=600,height=400)
      plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB (No.PC from HRC=",no.HRC,")",sep=""),ylab="g")
      points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
      points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
      points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
      legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
      dev.off()
      
    }
  }
}

for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.UKB in seq(0,100,5)){
      results2=results[which(results$income==income&results$location==location&results$PC_UKB==no.UKB),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB/income",income,"_location",location,"_no.UKB",no.UKB,"_HRC_varying.png",sep=""),width=600,height=400)
      plot(results2$PC_HRC,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from HRC (No.PC from UKB=",no.UKB,")",sep=""),ylab="g")
      points(results2$PC_HRC,results2$gc_all,pch=16,col="red")
      points(results2$PC_HRC,results2$gc_w.o.spo,pch=16,col="orange")
      points(results2$PC_HRC,results2$gc_spo,pch=16,col="red3")
      legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
      dev.off()
      
    }
  }
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
  year_check$plot_year[year_check$assigned_year==aa[idx]]=mean(bb)
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
  year_check_parents$plot_year[year_check_parents$assigned_year==aa[idx]]=mean(bb)
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
marriage_year_fat=marriage_year2[match(spo$IID_fat,marriage_year2$f.eid),]
marriage_year_mot=marriage_year2[match(spo$IID_mot,marriage_year2$f.eid),]
dir_name=paste("/pl/active/KellerLab/Yongkang/AM_analysis_results_PC40_trend/cov_include_year_sex_onlyEngland",sep="")    
dir_name2=paste("/pl/active/KellerLab/Yongkang/AM_analysis_results_PC40_trend/cov_include_year_sex_onlyEngland_Rdata",sep="")    
dir.create(dir_name)
dir.create(dir_name2)
spo_fat=NA;spo_mot=NA;bin_idx=NA



###############################Fit the model with 10-fold CV ####################################
library(data.table)
cor_inf=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_r2_cor_ver2.csv",header=T,stringsAsFactors=FALSE)
cor_inf=as.data.frame(cor_inf)


load("/pl/active/IBG/UKBiobank/PHENO/Keller/old_baskets/Keller-42049/ukb42049.RData")
race=data.frame(IID=bd[,1],race=bd[,11774])
british_IID=race$IID[which(race$race=="British")]
pheno2=pheno2[which(pheno2$IID%in%british_IID),]
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_phenofiles_no_scaled_addmore.csv",header=T,stringsAsFactors=FALSE)


bin_range_year=3
bin_range_std=5
bin_year=TRUE; bin_std=TRUE
phenonames=colnames(pheno2)
idx=13
phenoname=phenonames[idx]
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/",sep=""))%in%"chr_merged_results.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/chr_merged_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))%in%"chr_merged_results.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/chr_merged_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/",sep=""))%in%"chr_merged_results_ver2.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/chr_merged_results_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))%in%"chr_merged_results_ver2.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/chr_merged_results_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
#eigen_rare=fread("/pl/active/IBG/promero/rare-common-overlap/phenos_covar/rarePCs/pcs_thinnedPcs.txt",header=T,stringsAsFactors=FALSE) #Eigen vector from flash PCA with low frequency variables
#eigen_rare=as.data.frame(eigen_rare)

eigen_hrc=fread("/pl/active/KellerLab/Yongkang/HRC_pca/eigen_vec_from_hrc_to_UKB_top50.sscore",header=T,stringsAsFactors=FALSE)
eigen_hrc=as.data.frame(eigen_hrc)

eigen_vec=fread("/pl/active/IBG/UKBiobank/exclusions/relatedness/white/pca/projections.txt",header=T,stringsAsFactors=FALSE)  #Eigen vector from flash PCA 
eigen_vec=as.data.frame(eigen_vec)

results=data.frame(income=NA,location=NA,PC_UKB=NA,PC_HRC=NA,gt=NA,gc_all=NA,gc_w.o.spo=NA,gc_spo=NA)

for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.PC_HRC in seq(0,50,5)){
      for(no.PC_UKB in seq(0,100,5)){
        pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
        pheno2=as.data.frame(pheno2)
        covar=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_covarfile3.txt",header=T,stringsAsFactors=FALSE)
        covar=as.data.frame(covar)
        covar=covar[match(pheno2$IID,covar$id),]
        covar=covar[,-c(4:13)]
        marriage_year=fread("/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/Assortment_Project/Marriage_Years/UKB_SpousePairs_MarriageYears.txt",header=T,stringsAsFactors=FALSE)
        marriage_year=as.data.frame(marriage_year)
        marriage_year=marriage_year[match(pheno2$IID,marriage_year$f.eid),]
        marriage_year=marriage_year[,c(1,3,2,4)]
        spo=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/PGS/210517_NEW/spouse_dat_pheno_pgs_results.csv",header=T,stringsAsFactors=FALSE)
        spo=as.data.frame(spo)
        if(no.PC_HRC!=0&no.PC_UKB!=0){
          eigen_hrc2=eigen_hrc[match(pheno2$IID,eigen_hrc$IID),5:(4+no.PC_HRC)]
          eigen_vec2=eigen_vec[match(pheno2$IID,eigen_vec$IID),3:(2+no.PC_UKB)]
          eigen_vec2=cbind(eigen_hrc2,eigen_vec2)  
        }
        if(no.PC_HRC!=0&no.PC_UKB==0){
          eigen_hrc2=eigen_hrc[match(pheno2$IID,eigen_hrc$IID),5:(4+no.PC_HRC)]
          eigen_vec2=eigen_hrc2
        }
        if(no.PC_HRC==0&no.PC_UKB!=0){
          eigen_vec2=eigen_vec[match(pheno2$IID,eigen_vec$IID),3:(2+no.PC_UKB)]
        }
        
        if(income==0){
          covar=covar[,-23]
        }
        if(location==0){
          covar=covar[,-c(5:22)]
        }
        covar2=covar[,-1]
        if(no.PC_HRC!=0|no.PC_UKB!=0){
          covar2=cbind(covar2,eigen_vec2)
        }
        covar2$year_sex=covar2$Year_of_birth*covar2$Sex
        
        pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
        covar2=covar2[which(marriage_year$Participant.Marriage.Year<=1999),]
        marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]
        
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
        
        hpgs2=hpgs[match(pheno2$IID,hpgs$IID),]
        pgs=c(hpgs2[,2]+hpgs2[,3])
        
        samp=sample(1:length(pgs),length(pgs),replace=FALSE)
        
        #10-fold CV
        cv_term=as.integer(seq(0,length(pgs),length.out=11))
        for(idx in 1:10){
          pgs_train=pgs[-samp[(cv_term[idx]+1):cv_term[idx+1]]]
          covar_train=covar2[-samp[(cv_term[idx]+1):cv_term[idx+1]],]
          lmfit=lm(pgs_train~.,data=covar_train)
          covar_test=covar2[samp[(cv_term[idx]+1):cv_term[idx+1]],]
          fitted=predict(lmfit,newdata=covar_test)
          hpgs2[samp[(cv_term[idx]+1):cv_term[idx+1]],2]=hpgs2[samp[(cv_term[idx]+1):cv_term[idx+1]],2]-fitted/2
          hpgs2[samp[(cv_term[idx]+1):cv_term[idx+1]],3]=hpgs2[samp[(cv_term[idx]+1):cv_term[idx+1]],3]-fitted/2
          
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
        results2=data.frame(income=income,location=location,PC_UKB=no.PC_UKB,PC_HRC=no.PC_HRC,gt=gt,gc_all=gc_all,gc_w.o.spo=gc_w.o.spo,gc_spo=gc_spo)
        results=rbind(results,results2)
      }
    }
  }
}
results=results[-1,]
write.csv(results,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_compare_10fold_",phenoname,".csv",sep=""),row.names=F,quote=F)


###################Plot results (10-fold CV) #################################
results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_compare_10fold.csv",header=T,stringsAsFactors=FALSE)
results_nocv=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_compare.csv",header=T,stringsAsFactors=FALSE)

dir.create("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB/10foldCV_w.o.CV_combine")

for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.HRC in seq(0,50,5)){
      results2=results[which(results$income==income&results$location==location&results$PC_HRC==no.HRC),]
      results_nocv2=results_nocv[which(results_nocv$income==income&results_nocv$location==location&results_nocv$PC_HRC==no.HRC),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo,results_nocv2$gt,results_nocv2$gc_all,results_nocv2$gc_w.o.spo,results_nocv2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB/10foldCV_w.o.CV_combine/income",income,"_location",location,"_UKB_varying_no.HRC",no.HRC,".png",sep=""),width=600,height=400)
      plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB (No.PC from HRC=",no.HRC,")",sep=""),ylab="g",col="gray",type="b")
      points(results2$PC_UKB,results2$gc_all,pch=16,col="blue")
      lines(results2$PC_UKB,results2$gc_all,lty=2,col="blue")
      points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="skyblue")
      lines(results2$PC_UKB,results2$gc_w.o.spo,lty=2,col="skyblue")
      points(results2$PC_UKB,results2$gc_spo,pch=16,col="blue3")
      lines(results2$PC_UKB,results2$gc_spo,lty=2,col="blue3")
      points(results_nocv2$PC_UKB,results_nocv2$gt,pch=16)
      lines(results_nocv2$PC_UKB,results_nocv2$gt,lty=2)
      points(results_nocv2$PC_UKB,results_nocv2$gc_all,pch=16,col="red")
      lines(results_nocv2$PC_UKB,results_nocv2$gc_all,lty=2,col="red")
      points(results_nocv2$PC_UKB,results_nocv2$gc_w.o.spo,pch=16,col="orange")
      lines(results_nocv2$PC_UKB,results_nocv2$gc_w.o.spo,lty=2,col="orange")
      points(results_nocv2$PC_UKB,results_nocv2$gc_spo,pch=16,col="red3")
      lines(results_nocv2$PC_UKB,results_nocv2$gc_spo,col="red3",lty=2)
      legend("topright",pch=16,col=c("black","red","orange","red3","gray","blue","skyblue","blue3"),c("gt (No CV)","gc (all UKB,No CV)","gc (w.o. Spouses,No CV)","gc (only Spouses,No CV)","gt (10-fold CV)","gc (all UKB,10-fold CV)","gc (w.o. Spouses,10-fold CV)","gc (only Spouses,10-fold CV)"))
      dev.off()
      
    }
  }
}


for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.UKB in seq(0,100,5)){
      results2=results[which(results$income==income&results$location==location&results$PC_UKB==no.UKB),]
      results_nocv2=results_nocv[which(results_nocv$income==income&results_nocv$location==location&results_nocv$PC_UKB==no.UKB),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo,results_nocv2$gt,results_nocv2$gc_all,results_nocv2$gc_w.o.spo,results_nocv2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB/10foldCV_w.o.CV_combine/income",income,"_location",location,"_no.UKB",no.UKB,"_HRC_varying.png",sep=""),width=600,height=400)
      plot(results2$PC_HRC,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from HRC (No.PC from UKB=",no.UKB,")",sep=""),ylab="g",col="gray",type="b")
      points(results2$PC_HRC,results2$gc_all,pch=16,col="blue")
      lines(results2$PC_HRC,results2$gc_all,lty=2,col="blue")
      points(results2$PC_HRC,results2$gc_w.o.spo,pch=16,col="skyblue")
      lines(results2$PC_HRC,results2$gc_w.o.spo,lty=2,col="skyblue")
      points(results2$PC_HRC,results2$gc_spo,pch=16,col="blue3")
      lines(results2$PC_HRC,results2$gc_spo,lty=2,col="blue3")
      points(results_nocv2$PC_HRC,results_nocv2$gt,pch=16)
      lines(results_nocv2$PC_HRC,results_nocv2$gt,lty=2)
      points(results_nocv2$PC_HRC,results_nocv2$gc_all,pch=16,col="red")
      lines(results_nocv2$PC_HRC,results_nocv2$gc_all,lty=2,col="red")
      points(results_nocv2$PC_HRC,results_nocv2$gc_w.o.spo,pch=16,col="orange")
      lines(results_nocv2$PC_HRC,results_nocv2$gc_w.o.spo,lty=2,col="orange")
      points(results_nocv2$PC_HRC,results_nocv2$gc_spo,pch=16,col="red3")
      lines(results_nocv2$PC_HRC,results_nocv2$gc_spo,col="red3",lty=2)
      legend("topright",pch=16,col=c("black","red","orange","red3","gray","blue","skyblue","blue3"),c("gt (No CV)","gc (all UKB,No CV)","gc (w.o. Spouses,No CV)","gc (only Spouses,No CV)","gt (10-fold CV)","gc (all UKB,10-fold CV)","gc (w.o. Spouses,10-fold CV)","gc (only Spouses,10-fold CV)"))
      dev.off()
      
    }
  }
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
  year_check$plot_year[year_check$assigned_year==aa[idx]]=mean(bb)
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
  year_check_parents$plot_year[year_check_parents$assigned_year==aa[idx]]=mean(bb)
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
marriage_year_fat=marriage_year2[match(spo$IID_fat,marriage_year2$f.eid),]
marriage_year_mot=marriage_year2[match(spo$IID_mot,marriage_year2$f.eid),]
dir_name=paste("/pl/active/KellerLab/Yongkang/AM_analysis_results_PC40_trend/cov_include_year_sex_onlyEngland",sep="")    
dir_name2=paste("/pl/active/KellerLab/Yongkang/AM_analysis_results_PC40_trend/cov_include_year_sex_onlyEngland_Rdata",sep="")    
dir.create(dir_name)
dir.create(dir_name2)
spo_fat=NA;spo_mot=NA;bin_idx=NA

###################Plot results #################################
results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_compare_EA.csv",header=T,stringsAsFactors=FALSE)

dir.create("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB_EA")

for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.HRC in seq(0,50,5)){
      results2=results[which(results$income==income&results$location==location&results$PC_HRC==no.HRC),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB_EA/income",income,"_location",location,"_UKB_varying_no.HRC",no.HRC,".png",sep=""),width=600,height=400)
      plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB (No.PC from HRC=",no.HRC,")",sep=""),ylab="g",main="EA")
      points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
      points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
      points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
      legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
      dev.off()
      
    }
  }
}

for(income in c(0,1)){
  for(location in c(0,1)){
    for(no.UKB in seq(0,100,5)){
      results2=results[which(results$income==income&results$location==location&results$PC_UKB==no.UKB),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB_EA/income",income,"_location",location,"_no.UKB",no.UKB,"_HRC_varying.png",sep=""),width=600,height=400)
      plot(results2$PC_HRC,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from HRC (No.PC from UKB=",no.UKB,")",sep=""),ylab="g",main="EA")
      points(results2$PC_HRC,results2$gc_all,pch=16,col="red")
      points(results2$PC_HRC,results2$gc_w.o.spo,pch=16,col="orange")
      points(results2$PC_HRC,results2$gc_spo,pch=16,col="red3")
      legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
      dev.off()
      
    }
  }
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
  year_check$plot_year[year_check$assigned_year==aa[idx]]=mean(bb)
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
  year_check_parents$plot_year[year_check_parents$assigned_year==aa[idx]]=mean(bb)
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
marriage_year_fat=marriage_year2[match(spo$IID_fat,marriage_year2$f.eid),]
marriage_year_mot=marriage_year2[match(spo$IID_mot,marriage_year2$f.eid),]
dir_name=paste("/pl/active/KellerLab/Yongkang/AM_analysis_results_PC40_trend/cov_include_year_sex_onlyEngland",sep="")    
dir_name2=paste("/pl/active/KellerLab/Yongkang/AM_analysis_results_PC40_trend/cov_include_year_sex_onlyEngland_Rdata",sep="")    
dir.create(dir_name)
dir.create(dir_name2)
spo_fat=NA;spo_mot=NA;bin_idx=NA




#####################New Eigenvector each chromosome fitting for UKB samples #######
library(data.table)
cor_inf=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_r2_cor_ver2.csv",header=T,stringsAsFactors=FALSE)
cor_inf=as.data.frame(cor_inf)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno2=as.data.frame(pheno2)


load("/pl/active/IBG/UKBiobank/PHENO/Keller/old_baskets/Keller-42049/ukb42049.RData")
race=data.frame(IID=bd[,1],race=bd[,11774])
british_IID=race$IID[which(race$race=="British")]
pheno2=pheno2[which(pheno2$IID%in%british_IID),]
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_phenofiles_no_scaled_addmore.csv",header=T,stringsAsFactors=FALSE)

marriage_year=fread("/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/Assortment_Project/Marriage_Years/UKB_SpousePairs_MarriageYears.txt",header=T,stringsAsFactors=FALSE)
marriage_year=as.data.frame(marriage_year)
marriage_year=marriage_year[match(pheno2$IID,marriage_year$f.eid),]
marriage_year=marriage_year[,c(1,3,2,4)]
covar=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_covarfile3.txt",header=T,stringsAsFactors=FALSE)
covar=as.data.frame(covar)
covar=covar[match(pheno2$IID,covar$id),]
covar=covar[,-c(4:13)]
spo=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/PGS/210517_NEW/spouse_dat_pheno_pgs_results.csv",header=T,stringsAsFactors=FALSE)
spo=as.data.frame(spo)

bin_range_year=3
bin_range_std=5
bin_year=TRUE; bin_std=TRUE

phenoname="height"
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/",sep=""))%in%"chr_merged_results.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/chr_merged_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))%in%"chr_merged_results.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/chr_merged_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/",sep=""))%in%"chr_merged_results_ver2.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/chr_merged_results_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))%in%"chr_merged_results_ver2.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/chr_merged_results_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
#eigen_rare=fread("/pl/active/IBG/promero/rare-common-overlap/phenos_covar/rarePCs/pcs_thinnedPcs.txt",header=T,stringsAsFactors=FALSE) #Eigen vector from flash PCA with low frequency variables
#eigen_rare=as.data.frame(eigen_rare)

results=data.frame(income=NA,location=NA,PC_UKB=NA,gt=NA,gc_all=NA,gc_w.o.spo=NA,gc_spo=NA)
pheno2=pheno2[which(marriage_year$Participant.Marriage.Year<=1999),]
covar=covar[which(marriage_year$Participant.Marriage.Year<=1999),]
marriage_year=marriage_year[which(marriage_year$Participant.Marriage.Year<=1999),]

for(chr in 1:22){
  dat=fread(paste("/pl/active/KellerLab/Yongkang/UK_PC_each_chr/chr",chr,".eigenvec",sep=""),header=T,stringsAsFactors=FALSE)  #Eigen vector from flash PCA 
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
          hpgs_file=fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/individual_hpgs/chr",chr,".",idx,".csv",sep=""),header=T,stringsAsFactors=FALSE)
          hpgs_file=as.data.frame(hpgs_file)
          for(idx in 2:10){
            hpgs_file2=fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/individual_hpgs/chr",chr,".",idx,".csv",sep=""),header=T,stringsAsFactors=FALSE)
            hpgs_file2=as.data.frame(hpgs_file2)
            hpgs_file[,2:3]=hpgs_file[,2:3]+hpgs_file2[,2:3]
          }
          hpgs_file=hpgs_file[match(pheno2$IID,hpgs_file$IID),]
          covar3=covar2
          if(no.PC_UKB!=0){
            for(chr2 in (1:22)[-chr]){
            eigen_vec2=get(paste("chr",chr2,"_eigen",sep=""))[,3:(2+no.PC_UKB)]
            covar3=cbind(covar3,eigen_vec2)
            colnames(covar3)[(ncol(covar3)-no.PC_UKB+1):ncol(covar3)]=paste("PC_chr",chr2,"_",1:no.PC_UKB,sep="")
            }
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
        write.csv(results,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_compare_",phenoname,"_each_chr.csv",sep=""),row.names=F,quote=F)
        write.csv(R_squared,paste("/pl/active/KellerLab/Yongkang/HRC_pca/R_square_",phenoname,"_each_chr.csv",sep=""),row.names=F,quote=F)
      }
  }
}

results=results[-1,]
write.csv(results,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_compare_",phenoname,"_each_chr.csv",sep=""),row.names=F,quote=F)


#################PLotting results ###########
results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_compare_height_each_chr.csv",header=T,stringsAsFactors=FALSE)

dir.create("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB_EA_each_chr")
for(income in c(0,1)){
  for(location in c(0,1)){
      results2=results[which(results$income==income&results$location==location),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB_EA_each_chr/income",income,"_location",location,".png",sep=""),width=600,height=400)
      plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main="height")
      points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
      points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
      points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
      legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
      dev.off()
      
    }
  }

#####################New Eigenvector each chromosome fitting for UKB samples (Chr_merged_PC)#######
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
spo=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/PGS/210517_NEW/spouse_dat_pheno_pgs_results.csv",header=T,stringsAsFactors=FALSE)
spo=as.data.frame(spo)

bin_range_year=3
bin_range_std=5
bin_year=TRUE; bin_std=TRUE

phenoname="height"
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/",sep=""))%in%"chr_merged_results.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/chr_merged_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))%in%"chr_merged_results.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/chr_merged_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/",sep=""))%in%"chr_merged_results_ver2.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/whole_variant/chr_merged_results_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
if(sum(dir(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,sep=""))%in%"chr_merged_results_ver2.csv")==1){
  hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/chr_merged_results_ver2.csv",sep=""),header=T,stringsAsFactors=FALSE)
}
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
        hpgs_file=fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/individual_hpgs/chr",chr,".",idx,".csv",sep=""),header=T,stringsAsFactors=FALSE)
        hpgs_file=as.data.frame(hpgs_file)
        for(idx in 2:10){
          hpgs_file2=fread(paste("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/",phenoname,"/individual_hpgs/chr",chr,".",idx,".csv",sep=""),header=T,stringsAsFactors=FALSE)
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
      write.csv(results,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_compare_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)
      write.csv(R_squared,paste("/pl/active/KellerLab/Yongkang/HRC_pca/R_square_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)
    }
  }
}

results=results[-1,]
write.csv(results,paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_compare_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),row.names=F,quote=F)

#################PLotting results ###########
results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_compare_height_each_chr_PC_merged_chr.csv",header=T,stringsAsFactors=FALSE)
dir.create("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB_EA_each_chr_PC_merged_chr")
for(income in c(0,1)){
  for(location in c(0,1)){
    results2=results[which(results$income==income&results$location==location),]
    val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
    png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB_EA_each_chr_PC_merged_chr/income",income,"_location",location,".png",sep=""),width=600,height=400)
    plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main="height")
    points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
    points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
    points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
    legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
    dev.off()
    
  }
}

#################PLotting results ###########

phenonames=dir("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/new_PC_fitting/")
library(stringr)
phenonames=str_replace(phenonames,"R_square_","")
phenonames=str_replace(phenonames,"_each_chr_PC_merged_chr.csv","")
phenonames=str_replace(phenonames,"gc_gt_compare_","")
phenonames=unique(phenonames)
for(i in 1:length(phenonames)){
  phenoname=phenonames[i]
  results=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/new_PC_fitting/gc_gt_compare_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),header=T,stringsAsFactors=FALSE)
  for(income in c(0,1)){
    for(location in c(0,1)){
      results2=results[which(results$income==income&results$location==location),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/gc_gt_varying_no.UKB_EA_each_chr_PC_merged_chr/",phenoname,"_income",income,"_location",location,".png",sep=""),width=600,height=400)
      plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main=phenoname)
      points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
      points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
      points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
      legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
      dev.off()
      
    }
  }
  
}



phenonames=dir("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/new_PC_fitting_with_new_spouse_pair/")
dir.create("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/new_PC_fitting_with_new_spouse_pair/gc_gt_varying_no.UKB_EA_each_chr_PC_merged_chr")
library(stringr)
phenonames=str_replace(phenonames,"R_square_","")
phenonames=str_replace(phenonames,"_each_chr_PC_merged_chr.csv","")
phenonames=str_replace(phenonames,"gc_gt_compare_","")
phenonames=unique(phenonames)
phenonames=phenonames[phenonames!="gc_gt_varying_no.UKB_EA_each_chr_PC_merged_chr"]
for(i in 1:length(phenonames)){
  phenoname=phenonames[i]
  results=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/new_PC_fitting_with_new_spouse_pair//gc_gt_compare_",phenoname,"_each_chr_PC_merged_chr.csv",sep=""),header=T,stringsAsFactors=FALSE)
  for(income in c(0,1)){
    for(location in c(0,1)){
      results2=results[which(results$income==income&results$location==location),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/new_PC_fitting_with_new_spouse_pair/gc_gt_varying_no.UKB_EA_each_chr_PC_merged_chr/",phenoname,"_income",income,"_location",location,".png",sep=""),width=600,height=400)
      plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main=phenoname)
      points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
      points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
      points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
      legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
      dev.off()
      
    }
  }
  
}



results=read.csv(paste("/Users/yongkangkim/Downloads/gc_gt_compare_height_each_chr_PC_merged_chr.csv",sep=""),header=T,stringsAsFactors=FALSE)
for(income in c(0,1)){
  for(location in c(0,1)){
    results2=results[which(results$income==income&results$location==location),]
    val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
    png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/new_PC_fitting_with_new_spouse_pair/gc_gt_varying_no.UKB_EA_each_chr_PC_merged_chr/",phenoname,"_income",income,"_location",location,".png",sep=""),width=600,height=400)
    plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main=phenoname)
    points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
    points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
    points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
    legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
    dev.off()
    
  }
}

