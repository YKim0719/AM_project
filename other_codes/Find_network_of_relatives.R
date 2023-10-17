#####Read Pihat ##########
library(data.table)
dat=fread("/pl/active/KellerLab/Emmanuel/RunonfullUKB/deg123pihats.txt.genome",header=T,stringsAsFactors=FALSE)
dat=as.data.frame(dat)

###Initial round find relatives whose PI_hat >0.1##########
dat3=dat[dat$PI_HAT>0.1,]
ids2=c(dat3$IID1,dat3$IID2) #The list of all individuals who have at least one relatives in UKB
ta2=table(ids2) #Number of relatives for each individual
ta2=ta2[rev(order(ta2))] #Reordering the table with descending order
idx=as.numeric(names(ta2)[1]) #Pick individual who has the most number of relatives
dat4=dat3[which(dat3$IID1==idx|dat3$IID2==idx),] #Pairs which one of element is target individual
idx2=unique(c(dat4$IID1,dat4$IID2)) #The list of individuals who are related with picked individual
no.idx_former=length(idx)
  while(no.idx_former!=length(idx2)){
    no.idx_former=length(idx2)
    dat4=dat3[which(dat3$IID1%in%idx2|dat3$IID2%in%idx2),]
    idx2=unique(c(dat4$IID1,dat4$IID2))
  } #Find closed network including picked individual

no.relatives=data.frame(group=1,no.rel=length(idx2)) #Summary dataframe
group1_id=idx2 #List of individual for first closed network
dat3=dat3[-which(dat3$IID1%in%idx2|dat3$IID2%in%idx2),] #Exclude individuals in first closed network from Relative data frame 
while(nrow(dat3)!=0){
  ids2=c(dat3$IID1,dat3$IID2)
  ta2=table(ids2)
  ta2=ta2[rev(order(ta2))]
  
  idx=as.numeric(names(ta2)[1])
  dat4=dat3[which(dat3$IID1==idx|dat3$IID2==idx),]
  idx2=unique(c(dat4$IID1,dat4$IID2))
  no.idx_former=length(idx)
  while(no.idx_former!=length(idx2)){
    no.idx_former=length(idx2)
    dat4=dat3[which(dat3$IID1%in%idx2|dat3$IID2%in%idx2),]
    idx2=unique(c(dat4$IID1,dat4$IID2))
  }
  
  no.relatives2=data.frame(group=nrow(no.relatives)+1,no.rel=length(idx2))
  no.relatives=rbind(no.relatives,no.relatives2)
  assign(paste("group",nrow(no.relatives),"_id",sep=""),idx2)
  dat3=dat3[-which(dat3$IID1%in%idx2|dat3$IID2%in%idx2),]
}
dir.create("/pl/active/KellerLab/Yongkang/list_of_relative_network/")
save(file="/pl/active/KellerLab/Yongkang/list_of_relative_network/group_ids_collection_pihat0.1.RData",list=paste("group",nrow(no.relatives),"_id",sep=""))
write.csv(no.relatives,"/pl/active/KellerLab/Yongkang/list_of_relative_network/report_no.relatives_each_group_pihat0.1.csv",row.names=F,quote=F)


###Initial round finding relatives whose PI_hat >0.05##########
dat3=dat[dat$PI_HAT>0.05,]
ids2=c(dat3$IID1,dat3$IID2)
ta2=table(ids2)
ta2=ta2[rev(order(ta2))]

idx=as.numeric(names(ta2)[1])
dat4=dat3[which(dat3$IID1==idx|dat3$IID2==idx),]
idx2=unique(c(dat4$IID1,dat4$IID2))
no.idx_former=length(idx)
while(no.idx_former!=length(idx2)){
  no.idx_former=length(idx2)
  dat4=dat3[which(dat3$IID1%in%idx2|dat3$IID2%in%idx2),]
  idx2=unique(c(dat4$IID1,dat4$IID2))
}

no.relatives=data.frame(group=1,no.rel=length(idx2))
group1_id=idx2
dat3=dat3[-which(dat3$IID1%in%idx2|dat3$IID2%in%idx2),]
while(nrow(dat3)!=0){
  ids2=c(dat3$IID1,dat3$IID2)
  ta2=table(ids2)
  ta2=ta2[rev(order(ta2))]
  
  idx=as.numeric(names(ta2)[1])
  dat4=dat3[which(dat3$IID1==idx|dat3$IID2==idx),]
  idx2=unique(c(dat4$IID1,dat4$IID2))
  no.idx_former=length(idx)
  while(no.idx_former!=length(idx2)){
    no.idx_former=length(idx2)
    dat4=dat3[which(dat3$IID1%in%idx2|dat3$IID2%in%idx2),]
    idx2=unique(c(dat4$IID1,dat4$IID2))
  }
  
  no.relatives2=data.frame(group=nrow(no.relatives)+1,no.rel=length(idx2))
  no.relatives=rbind(no.relatives,no.relatives2)
  assign(paste("group",nrow(no.relatives),"_id",sep=""),idx2)
  dat3=dat3[-which(dat3$IID1%in%idx2|dat3$IID2%in%idx2),]
}
save(file="/pl/active/KellerLab/Yongkang/list_of_relative_network/group_ids_collection_pihat0.05.RData",list=paste("group",nrow(no.relatives),"_id",sep=""))
write.csv(no.relatives,"/pl/active/KellerLab/Yongkang/list_of_relative_network/report_no.relatives_each_group_pihat0.05.csv",row.names=F,quote=F)
