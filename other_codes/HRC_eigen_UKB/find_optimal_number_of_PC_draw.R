library(stringr)
aa=dir("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/results_internal_analysis/optimal_hpgs/")
aa=str_replace(aa,"gc_gt_compare_","")
aa=str_replace(aa,"_each_chr_PC_merged_chr.csv","")
aa=str_replace(aa,"R_square_","")
aa=unique(aa)

idx=1
dat=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/results_internal_analysis/optimal_hpgs/gc_gt_compare_",aa[idx],"_each_chr_PC_merged_chr.csv",sep=""),header=T,stringsAsFactors=FALSE)
dir.create("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/results_internal_analysis/figures")

for(idx in 1:length(aa)){
  results=read.csv(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/results_internal_analysis/optimal_hpgs/gc_gt_compare_",aa[idx],"_each_chr_PC_merged_chr.csv",sep=""),header=T,stringsAsFactors=FALSE)
  phenoname=aa[idx]
  dir.create(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/results_internal_analysis/figures/",phenoname,sep=""))
  for(income in c(0,1)){
    for(location in c(0,1)){
      results2=results[which(results$income==income&results$location==location),]
      val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
      png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/results_internal_analysis/figures/",phenoname,"/income",income,"_location",location,".png",sep=""),width=600,height=400)
      plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main=phenoname)
      points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
      points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
      points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
      legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
      dev.off()
      
    }
  }
  
}
#Use array to perform bunch of jobs 


#######################Check ###############
results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/PC_check/external_gc_gt_compare_height_each_chr_PC_merged_chr.csv",header=T,stringsAsFactors=FALSE)
income=1;location=1
results2=results[which(results$income==income&results$location==location),]
val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
phenoname="height"
png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/PC_check/external_gc_gt_exclude_chr_dat.png",sep=""),width=600,height=400)
plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main=phenoname)
points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
dev.off()

results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/PC_check/external_gc_gt_compare_height_each_chr_PC_merged_chr_use_chr_PC.csv",header=T,stringsAsFactors=FALSE)
income=1;location=1
results2=results[which(results$income==income&results$location==location),]
val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
phenoname="height"
png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/PC_check/external_gc_gt_include_chr_dat.png",sep=""),width=600,height=400)
plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main=phenoname)
points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
dev.off()


results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/PC_check/internal_gc_gt_compare_height_each_chr_PC_merged_chr.csv",header=T,stringsAsFactors=FALSE)
income=1;location=1
results2=results[which(results$income==income&results$location==location),]
val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
phenoname="height"
png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/PC_check/internal_gc_gt_exclude_chr_dat.png",sep=""),width=600,height=400)
plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main=phenoname)
points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
legend("right",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
dev.off()

results=read.csv("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/PC_check/internal_gc_gt_compare_height_each_chr_PC_merged_chr_use_chr_PC.csv",header=T,stringsAsFactors=FALSE)
income=1;location=1
results2=results[which(results$income==income&results$location==location),]
val=c(results2$gt,results2$gc_all,results2$gc_w.o.spo,results2$gc_spo)
phenoname="height"
png(paste("/Users/yongkangkim/OneDrive - UCB-O365/Matthew c keller/sibling_VT_SEM/analysis_2023/PC_check/internal_gc_gt_include_chr_dat.png",sep=""),width=600,height=400)
plot(results2$PC_UKB,results2$gt,pch=16,ylim=c(min(val),max(val)),xlab=paste("No.PC from UKB)",sep=""),ylab="g",main=phenoname)
points(results2$PC_UKB,results2$gc_all,pch=16,col="red")
points(results2$PC_UKB,results2$gc_w.o.spo,pch=16,col="orange")
points(results2$PC_UKB,results2$gc_spo,pch=16,col="red3")
legend("topright",pch=16,col=c("black","red","orange","red3"),c("gt","gc (all UKB)","gc (w.o. Spouses)","gc (only Spouses"))
dev.off()

