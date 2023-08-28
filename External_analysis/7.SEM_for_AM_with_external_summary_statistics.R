args <- commandArgs(trailingOnly =TRUE)
idx=as.numeric(args[1])

###Cleaning Codes to make reproducible results

##singularity run /projects/lessem/singularity/openmx.sif Rscript /pl/active/KellerLab/Yongkang/AM_source_code/find_SEM_results_include_year_sex.R $idx

#####Generate the new external summary statistic bassed results ####
library(data.table)
source("/pl/active/KellerLab/Yongkang/AM_source_code/source_code_for_SEM_ver2.R")
source("/pl/active/KellerLab/Yongkang/AM_source_code/source_code_for_various_hpgs_handle_ver2.R")
cor_inf=fread("/rc_scratch/yoki5348/UK_PC_each_chr/gc_gt_internal/R_sqaured_internal.csv",header=T,stringsAsFactors=FALSE)
cor_inf=as.data.frame(cor_inf)
#pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno2=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_british.txt",header=T,stringsAsFactors=FALSE)
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

hpgs=read.csv(paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/hpgs_collection/",phenoname,"/hpgs_pc25.csv",sep=""),header=T,stringsAsFactors=FALSE)

hpgs2=hpgs[match(pheno2$IID,hpgs$IID),]

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
dir_name=paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/SEM_analysis_results/PC25",sep="")    
dir_name2=paste("/pl/active/KellerLab/Yongkang/HRC_pca/gc_gt_new_spouse/SEM_analysis_results/PC25_Rdata",sep="")    
dir.create(dir_name)
dir.create(dir_name2)
spo_fat=NA;spo_mot=NA;bin_idx=NA

pheno_fat=pheno_dat[match(spo$IID_fat,pheno_dat$IID),2]
pheno_mot=pheno_dat[match(spo$IID_mot,pheno_dat$IID),2]

if(phenoname%in%c("schizophrenia","bipolar","pgs_major_depression_non_UKB","anxiety_cc","ocd","PTSD","alzheimer")){
  g_result=try(perform_SEM_gc_gt_same_with_conti_w.o.y(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  global_g_result=try(perform_SEM_global_with_conti_w.o.y(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  trend_result=try(perform_SEM_linear_trend_with_conti_w.o.y(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  nonlinear_trend_result=try(perform_SEM_nonlinear_trend_with_conti_w.o.y(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  cor_ref=cor_inf$cor_ref[cor_inf$pheno==phenoname]
  #Omega2mu_gc_gt=mxEval(4*gt1/se_ref-cor_ref,global_g_result)
  #Omega2mu_gc_gt.se=mxSE(4*gt1/se_ref-cor_ref,global_g_result)
  gtest_1967=mxEval(2*gc11-gc1,nonlinear_trend_result)
  gtest_1967.se=mxSE(2*gc11-gc1,nonlinear_trend_result)
  gcest_1937=mxEval((gt1+gc11)/2,nonlinear_trend_result)
  gcest_1937.se= mxSE((gt1+gc11)/2,nonlinear_trend_result)
  rhoest_1967=mxEval(2*gcest11-gcest1,nonlinear_trend_result)
  rhoest_1967.se=mxSE(2*gcest11-gcest1,nonlinear_trend_result)
  rhoest_1937=mxEval((gtest1+gcest11)/2,nonlinear_trend_result)
  rhoest_1937.se= mxSE((gtest1+gcest11)/2,nonlinear_trend_result)
}else{
  g_result=try(perform_SEM_gc_gt_same_with_conti(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  global_g_result=try(perform_SEM_global_with_conti(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  trend_result=try(perform_SEM_linear_trend_with_conti(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
  nonlinear_trend_result=try(perform_SEM_nonlinear_trend_with_conti(phenoname,hpgs2,covar2,pheno_dat,cor_inf,dir_name))
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
write.csv(bb,paste(dir_name,"/",phenoname,"/lrt_test_results.csv",sep=""),row.names=F,quote=F)

est_trend=data.frame(phenoname,gtest_1967,
                     gtest_1967.se,
                     gcest_1937,
                     gcest_1937.se,
                     rhoest_1967,
                     rhoest_1967.se,
                     rhoest_1937,
                     rhoest_1937.se
                     #Omega2mu_gc_gt,
                     #Omega2mu_gc_gt.se
)
write.csv(est_trend,paste(dir_name,"/",phenoname,"/est_trend.csv",sep=""),row.names=F,quote=F)
dir.create(paste(dir_name2,"/",phenoname,sep=""))
save.image(file=paste(dir_name2,"/",phenoname,"/Rimage.RData",sep=""))

dat1=read.csv(paste(dir_name,"/",phenoname,"/diseq_AM_assumption_gc_gt_nonlinear.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat2=read.csv(paste(dir_name,"/",phenoname,"/diseq_AM_assumption_gc_gt_linear_increasing.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat3=read.csv(paste(dir_name,"/",phenoname,"/diseq_AM_assumption_gc_gt_global.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat4=read.csv(paste(dir_name,"/",phenoname,"/diseq_AM_assumption_gc_gt_same.csv",sep=""),header=T,stringsAsFactors=FALSE)
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


png(paste(dir_name,"/",phenoname,"/rhoplot_fiex_axis.png",sep=""),width=1500,height=800,pointsize=16)
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




#####Generate the new results for internal analysis



########Generate simulation results with PCs without target chromosomes

##########Check relationship between PCs and geographical information & income