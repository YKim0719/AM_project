library("plotrix")
aa=dir("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/")
phenonames=c("height","BMI","WHR","body_fat","EA","nervous","bdd","life_mdd","AgeSmk","CigDay_include_never_smoke","SmkCes","SmkInit","DrnkWk","diastolic_blood_pressure","systolic_blood_pressure","tg","hdl","ldl","T2D")
phenonames_init=c("Height","BMI","WHR","Fat","EA","NERV","BDD","MDD","AgeSmk","CigDay","SmkCes","SmkInit","DrnkWk","DBP","SBP","TG","HDL","LDL","T2D")
dat1=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/TG/diseq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat1[3,5]
phenonames_corinf=c("height","BMI","WHR","body_fat","EA","pgs_major_depression_non_UKB","anxiety_cc" ,"BDD" ,
                    "AgeOfInitiation_gscan","CigarettesPerDay_gscan" ,"SmokingCessation_gscan" ,"SmokingInitiation_gscan" ,
                    "DrinksPerWeek_gscan","DBP","SBP","TG","HDL","LDL","T2D")
cor_inf=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/UKB_pheno_r2_cor.csv",header=T,stringsAsFactors=FALSE)
cor_inf2=data.frame(pheno="BDD",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.1538898)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="body_fat",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.139388)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="T2D",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.2044781)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="TG",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.15035045)
cor_inf=rbind(cor_inf,cor_inf2)

cor_inf2=cor_inf[match(phenonames_corinf,cor_inf[,1]),]

cor_inf2$cor_ref[c(9,14,15)]=c(0.10100000,0.18422089,0.07541558)
idx=1
phenoname=phenonames[idx]
pheno_used=phenonames_init[idx]
dat1=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/diseq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat2=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/eq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
if(sum(dat1$name=="gtest")==1){
  gtest=dat1[dat1$name=="gtest",5]
  gtest_lower.C.I=dat1[dat1$name=="gtest",7]
  gtest_upper.C.I=dat1[dat1$name=="gtest",8]
  gcest=dat1[dat1$name=="gcest",5]
  gcest_lower.C.I=dat1[dat1$name=="gcest",7]
  gcest_upper.C.I=dat1[dat1$name=="gcest",8]
  gest=dat2[dat2$name=="rhogest",5]
  gest_lower.C.I=dat2[dat2$name=="rhogest",7]
  gest_upper.C.I=dat2[dat2$name=="rhogest",8]
  
  results=data.frame(pheno=pheno_used,rho_gc=gcest,rho_gc.li=gcest_lower.C.I,rho_gc.ui=gcest_upper.C.I,
                     rho_gt=gtest,rho_gt.li=gtest_lower.C.I,rho_gt.ui=gtest_upper.C.I,
                     rho_g=gest,rho_g.li=gest_lower.C.I,rho_g.ui=gest_upper.C.I
  )
  
}else{
  gtest_ref=dat1[dat1$name=="gtest_ref",5]
  gtest_ref_lower.C.I=dat1[dat1$name=="gtest_ref",7]
  gtest_ref_lower.C.I=dat1[dat1$name=="gtest_ref",8]
  gcest_ref=dat1[dat1$name=="gces_reft",5]
  gcest_ref_lower.C.I=dat1[dat1$name=="gcest_ref",7]
  gcest_ref_lower.C.I=dat1[dat1$name=="gcest_ref",8]
  gest_ref=dat2[dat2$name=="rhogest_ref",5]
  gest_ref_lower.C.I=dat2[dat2$name=="rhogest_ref",7]
  gest_ref_upper.C.I=dat2[dat2$name=="rhogest_ref",8]
  results=data.frame(pheno=pheno_used,rho_gc=gcest_ref,rho_gc.li=gcest_ref_lower.C.I,rho_gc.ui=gcest_ref_upper.C.I,
                     rho_gt=gtest_ref,rho_gt.li=gtest_ref_lower.C.I,rho_gt.ui=gtest_ref_upper.C.I,
                     rho_g=gest_ref,rho_g.li=gest_ref_lower.C.I,rho_g.ui=gest_ref_upper.C.I
  )
  
}

for(idx in 2:length(phenonames)){
  phenoname=phenonames[idx]
  pheno_used=phenonames_init[idx]
  dat1=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/diseq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
  dat2=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/eq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
  if(sum(dat1$name=="gtest")==1){
    gtest=dat1[dat1$name=="gtest",5]
    gtest_lower.C.I=dat1[dat1$name=="gtest",7]
    gtest_upper.C.I=dat1[dat1$name=="gtest",8]
    gcest=dat1[dat1$name=="gcest",5]
    gcest_lower.C.I=dat1[dat1$name=="gcest",7]
    gcest_upper.C.I=dat1[dat1$name=="gcest",8]
    gest=dat2[dat2$name=="rhogest",5]
    gest_lower.C.I=dat2[dat2$name=="rhogest",7]
    gest_upper.C.I=dat2[dat2$name=="rhogest",8]
    
    results2=data.frame(pheno=pheno_used,rho_gc=gcest,rho_gc.li=gcest_lower.C.I,rho_gc.ui=gcest_upper.C.I,
                       rho_gt=gtest,rho_gt.li=gtest_lower.C.I,rho_gt.ui=gtest_upper.C.I,
                       rho_g=gest,rho_g.li=gest_lower.C.I,rho_g.ui=gest_upper.C.I
    )
    
  }else{
    gtest_ref=dat1[dat1$name=="gtest_ref",5]
    gtest_ref_lower.C.I=dat1[dat1$name=="gtest_ref",7]
    gtest_ref_lower.C.I=dat1[dat1$name=="gtest_ref",8]
    gcest_ref=dat1[dat1$name=="gces_reft",5]
    gcest_ref_lower.C.I=dat1[dat1$name=="gcest_ref",7]
    gcest_ref_lower.C.I=dat1[dat1$name=="gcest_ref",8]
    gest_ref=dat2[dat2$name=="rhogest_ref",5]
    gest_ref_lower.C.I=dat2[dat2$name=="rhogest_ref",7]
    gest_ref_upper.C.I=dat2[dat2$name=="rhogest_ref",8]
    results2=data.frame(pheno=pheno_used,rho_gc=gcest_ref,rho_gc.li=gcest_ref_lower.C.I,rho_gc.ui=gcest_ref_upper.C.I,
                       rho_gt=gtest_ref,rho_gt.li=gtest_ref_lower.C.I,rho_gt.ui=gtest_ref_upper.C.I,
                       rho_g=gest_ref,rho_g.li=gest_ref_lower.C.I,rho_g.ui=gest_ref_upper.C.I
    )
    
  }
  
  results=rbind(results,results2)
}


results_plot=data.frame(x=rep(1:nrow(results),3)+rep(c(0,0.2,0.4),each=nrow(results)),y=c(results$rho_gc,results$rho_gt,results$rho_g),li=c(results$rho_gc.li,results$rho_gt.li,results$rho_g.li),ui=c(results$rho_gc.ui,results$rho_gt.ui,results$rho_g.ui))
plotCI(x = results_plot$x,y =results_plot$y,li = results_plot$li,ui = results_plot$ui,axes=FALSE,scol=rep(c("skyblue","red","purple"),each=nrow(results)),
       xlab="",ylab="Estimates",ylim=c(-0.1,2))#scol sfrac=1
for(id in 1:nrow(cor_inf2)){
  lines(c(id,id+0.4),c(cor_inf2$cor_ref[id],cor_inf2$cor_ref[id]),lwd=2,lty=2,col="green")
}

abline(a=0,b=0,lty=2)
axis(1,phenonames_init,at=c(1:nrow(results))+0.2)
axis(2)
box()


######################Selected phenotype####################
phenonames=c("height","BMI","WHR","EA","bdd","SmkCes","SmkInit","DrnkWk","diastolic_blood_pressure","systolic_blood_pressure","tg","hdl","ldl","T2D")
phenonames_init=c("Height","BMI","WHR","EA","BDD","SmkCes","SmkInit","DrnkWk","DBP","SBP","TG","HDL","LDL","T2D")
dat1=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/TG/diseq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat1[3,5]
phenonames_corinf=c("height","BMI","WHR","EA" ,"BDD" ,
                    "SmokingCessation_gscan" ,"SmokingInitiation_gscan" ,
                    "DrinksPerWeek_gscan","DBP","SBP","TG","HDL","LDL","T2D")
cor_inf=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/UKB_pheno_r2_cor.csv",header=T,stringsAsFactors=FALSE)
cor_inf2=data.frame(pheno="BDD",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.1538898)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="body_fat",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.139388)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="T2D",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.2044781)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="TG",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.15035045)
cor_inf=rbind(cor_inf,cor_inf2)

cor_inf2=cor_inf[match(phenonames_corinf,cor_inf[,1]),]

cor_inf2$cor_ref[c(9,10)]=c(0.18422089,0.07541558)
idx=1
phenoname=phenonames[idx]
pheno_used=phenonames_init[idx]
dat1=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/diseq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
dat2=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/eq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
if(sum(dat1$name=="gtest")==1){
  gtest=dat1[dat1$name=="gtest",5]
  gtest_lower.C.I=dat1[dat1$name=="gtest",7]
  gtest_upper.C.I=dat1[dat1$name=="gtest",8]
  gcest=dat1[dat1$name=="gcest",5]
  gcest_lower.C.I=dat1[dat1$name=="gcest",7]
  gcest_upper.C.I=dat1[dat1$name=="gcest",8]
  gest=dat2[dat2$name=="rhogest",5]
  gest_lower.C.I=dat2[dat2$name=="rhogest",7]
  gest_upper.C.I=dat2[dat2$name=="rhogest",8]
  
  results=data.frame(pheno=pheno_used,rho_gc=gcest,rho_gc.li=gcest_lower.C.I,rho_gc.ui=gcest_upper.C.I,
                     rho_gt=gtest,rho_gt.li=gtest_lower.C.I,rho_gt.ui=gtest_upper.C.I,
                     rho_g=gest,rho_g.li=gest_lower.C.I,rho_g.ui=gest_upper.C.I
  )
  
}else{
  gtest_ref=dat1[dat1$name=="gtest_ref",5]
  gtest_ref_lower.C.I=dat1[dat1$name=="gtest_ref",7]
  gtest_ref_lower.C.I=dat1[dat1$name=="gtest_ref",8]
  gcest_ref=dat1[dat1$name=="gces_reft",5]
  gcest_ref_lower.C.I=dat1[dat1$name=="gcest_ref",7]
  gcest_ref_lower.C.I=dat1[dat1$name=="gcest_ref",8]
  gest_ref=dat2[dat2$name=="rhogest_ref",5]
  gest_ref_lower.C.I=dat2[dat2$name=="rhogest_ref",7]
  gest_ref_upper.C.I=dat2[dat2$name=="rhogest_ref",8]
  results=data.frame(pheno=pheno_used,rho_gc=gcest_ref,rho_gc.li=gcest_ref_lower.C.I,rho_gc.ui=gcest_ref_upper.C.I,
                     rho_gt=gtest_ref,rho_gt.li=gtest_ref_lower.C.I,rho_gt.ui=gtest_ref_upper.C.I,
                     rho_g=gest_ref,rho_g.li=gest_ref_lower.C.I,rho_g.ui=gest_ref_upper.C.I
  )
  
}

for(idx in 2:length(phenonames)){
  phenoname=phenonames[idx]
  pheno_used=phenonames_init[idx]
  dat1=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/diseq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
  dat2=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/eq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
  if(sum(dat1$name=="gtest")==1){
    gtest=dat1[dat1$name=="gtest",5]
    gtest_lower.C.I=dat1[dat1$name=="gtest",7]
    gtest_upper.C.I=dat1[dat1$name=="gtest",8]
    gcest=dat1[dat1$name=="gcest",5]
    gcest_lower.C.I=dat1[dat1$name=="gcest",7]
    gcest_upper.C.I=dat1[dat1$name=="gcest",8]
    gest=dat2[dat2$name=="rhogest",5]
    gest_lower.C.I=dat2[dat2$name=="rhogest",7]
    gest_upper.C.I=dat2[dat2$name=="rhogest",8]
    
    results2=data.frame(pheno=pheno_used,rho_gc=gcest,rho_gc.li=gcest_lower.C.I,rho_gc.ui=gcest_upper.C.I,
                        rho_gt=gtest,rho_gt.li=gtest_lower.C.I,rho_gt.ui=gtest_upper.C.I,
                        rho_g=gest,rho_g.li=gest_lower.C.I,rho_g.ui=gest_upper.C.I
    )
    
  }else{
    gtest_ref=dat1[dat1$name=="gtest_ref",5]
    gtest_ref_lower.C.I=dat1[dat1$name=="gtest_ref",7]
    gtest_ref_lower.C.I=dat1[dat1$name=="gtest_ref",8]
    gcest_ref=dat1[dat1$name=="gces_reft",5]
    gcest_ref_lower.C.I=dat1[dat1$name=="gcest_ref",7]
    gcest_ref_lower.C.I=dat1[dat1$name=="gcest_ref",8]
    gest_ref=dat2[dat2$name=="rhogest_ref",5]
    gest_ref_lower.C.I=dat2[dat2$name=="rhogest_ref",7]
    gest_ref_upper.C.I=dat2[dat2$name=="rhogest_ref",8]
    results2=data.frame(pheno=pheno_used,rho_gc=gcest_ref,rho_gc.li=gcest_ref_lower.C.I,rho_gc.ui=gcest_ref_upper.C.I,
                        rho_gt=gtest_ref,rho_gt.li=gtest_ref_lower.C.I,rho_gt.ui=gtest_ref_upper.C.I,
                        rho_g=gest_ref,rho_g.li=gest_ref_lower.C.I,rho_g.ui=gest_ref_upper.C.I
    )
    
  }
  
  results=rbind(results,results2)
}


results_plot=data.frame(x=rep(1:nrow(results),3)+rep(c(0,0.2,0.4),each=nrow(results)),y=c(results$rho_gc,results$rho_gt,results$rho_g),li=c(results$rho_gc.li,results$rho_gt.li,results$rho_g.li),ui=c(results$rho_gc.ui,results$rho_gt.ui,results$rho_g.ui))
plotCI(x = results_plot$x,y =results_plot$y,li = results_plot$li,ui = results_plot$ui,axes=FALSE,scol=rep(c("skyblue","red","purple"),each=nrow(results)),
       xlab="",ylab="Estimates",ylim=c(-0.1,2))#scol sfrac=1
for(id in 1:nrow(cor_inf2)){
  lines(c(id,id+0.4),c(cor_inf2$cor_ref[id],cor_inf2$cor_ref[id]),lwd=2,lty=2,col="green")
}

abline(a=0,b=0,lty=2)
axis(1,phenonames_init,at=c(1:nrow(results))+0.2)
axis(2)
box()


############LRT results #######
phenonames=c("height","BMI","WHR","body_fat","EA","nervous","bdd","life_mdd","AgeSmk","CigDay_include_never_smoke","SmkCes","SmkInit","DrnkWk","diastolic_blood_pressure","systolic_blood_pressure","tg","hdl","ldl","T2D")
phenonames_init=c("Height","BMI","WHR","Fat","EA","NERV","BDD","MDD","AgeSmk","CigDay","SmkCes","SmkInit","DrnkWk","DBP","SBP","TG","HDL","LDL","T2D")
phenonames_corinf=c("height","BMI","WHR","body_fat","EA","pgs_major_depression_non_UKB","anxiety_cc" ,"BDD" ,
                    "AgeOfInitiation_gscan","CigarettesPerDay_gscan" ,"SmokingCessation_gscan" ,"SmokingInitiation_gscan" ,
                    "DrinksPerWeek_gscan","DBP","SBP","TG","HDL","LDL","T2D")
cor_inf=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/results/UKB_pheno_r2_cor.csv",header=T,stringsAsFactors=FALSE)
cor_inf2=data.frame(pheno="BDD",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.1538898)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="body_fat",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.139388)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="T2D",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.2044781)
cor_inf=rbind(cor_inf,cor_inf2)
cor_inf2=data.frame(pheno="TG",R2_UKB=NA,R2_ref=NA,cor_UKB=NA,cor_ref=0.15035045)
cor_inf=rbind(cor_inf,cor_inf2)

cor_inf2=cor_inf[match(phenonames_corinf,cor_inf[,1]),]

cor_inf2$cor_ref[c(9,14,15)]=c(0.10100000,0.18422089,0.07541558)

lrt=c()
phenotypic_homogamy=c()
signs=c()
for(idx in 1:length(phenonames)){
  phenoname=phenonames[idx]
  pheno_used=phenonames_init[idx]
  dat1=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/lrt_test_results.csv",sep=""),header=T,stringsAsFactors=FALSE)
  lrt=c(lrt,dat1$p[2])
  if(dat1$p[2]<=0.05){
    dat2=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/diseq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
    if(sum(dat2$name=="gtest")==1){
      phenotypic_homogamy=c(phenotypic_homogamy,2*(1-pnorm(abs(dat2$Estimate[dat2$name=="gtest"]-cor_inf2$cor_ref[idx])/dat2$Std.Error[dat2$name=="gtest"])))
      signs=c(signs,dat2$Estimate[dat2$name=="gtest"]-cor_inf2$cor_ref[idx])
    }else{
      phenotypic_homogamy=c(phenotypic_homogamy,2*(1-pnorm(abs(dat2$Estimate[dat2$name=="gtest_ref"]-cor_inf2$cor_ref[idx])/dat2$Std.Error[dat2$name=="gtest_ref"])))
      signs=c(signs,dat2$Estimate[dat2$name=="gtest_ref"]-cor_inf2$cor_ref[idx])
      
    }
  }else{
    dat2=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/eq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
    if(sum(dat2$name=="rhogest")==1){
      phenotypic_homogamy=c(phenotypic_homogamy,2*(1-pnorm(abs(dat2$Estimate[dat2$name=="rhogest"]-cor_inf2$cor_ref[idx])/dat2$Std.Error[dat2$name=="rhogest"])))
      signs=c(signs,dat2$Estimate[dat2$name=="rhogest"]-cor_inf2$cor_ref[idx])
      
    }else{
      phenotypic_homogamy=c(phenotypic_homogamy,2*(1-pnorm(abs(dat2$Estimate[dat2$name=="rhogest_ref"]-cor_inf2$cor_ref[idx])/dat2$Std.Error[dat2$name=="rhogest_ref"])))
      signs=c(signs,dat2$Estimate[dat2$name=="rhogest_ref"]-cor_inf2$cor_ref[idx])
      
    }
    
  }
}
results=data.frame(pheno=phenonames,gc_gt_equal=lrt,homogamy=phenotypic_homogamy,effects=signs)
write.csv(results,"/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/inference_results.csv",row.names=F,quote=F)
  dat2=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/diseq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
  dat3=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/SEM_results_no_trend/",phenoname,"/eq_AM_assumption.csv",sep=""),header=T,stringsAsFactors=FALSE)
  
