##### Read Pihat ##########
library(data.table)
dat=fread("/pl/active/KellerLab/Emmanuel/RunonfullUKB/deg123pihats.txt.genome",header=T,stringsAsFactors=FALSE)
dat=as.data.frame(dat)
spo=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/PGS/210517_NEW/spouse_dat_pheno_pgs_results.csv",header=T,stringsAsFactors=FALSE)
spo=as.data.frame(spo)
##Find samples who have at least one relative (pi hat >0.05) or spouse among UKB paritcipants ##########

dat2=dat[dat$PI_HAT>0.05,]
dat2=data.frame(IID1=dat2$IID1,IID2=dat2$IID2)
spo2=data.frame(IID1=spo$IID_fat,IID2=spo$IID_mot)
dat2=rbind(dat2,spo2)
#dat2=dat2[which((dat2$IID1%in%spo$IID_fat& dat2$IID2%in%spo$IID_mot)|(dat2$IID1%in%spo$IID_mot&dat2$IID2%in%spo$IID_fat)),]
idx=1
aa=c(dat2[,1],dat2[,2])
aa2=unique(aa)
write.table(aa2,"/projects/yoki5348/test_samp_ID_all_spo.txt",row.names=F,col.names=F,quote=F)
library(data.table)
pheno=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno=as.data.frame(pheno)
bb=pheno$IID[-which(pheno$IID%in%aa2)]
write.table(bb,"/projects/yoki5348/train_samp_ID_w.o.spo.txt",row.names=F,col.names=F,quote=F)

train_ids=fread("/projects/yoki5348/train_samp_ID_w.o.spo.txt",header=F,stringsAsFactors=FALSE)
train_ids=as.data.frame(train_ids)[,1]
pheno2=pheno[match(train_ids,pheno$IID),]
write.table(pheno2,"/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_train.txt",row.names=F,col.names=T,quote=F,sep="\t")

LA.CA <- read.csv("UKB.AMC.current_address.UK_local_authorities.csv", header=T, stringsAsFactors = FALSE)
MSOA.CA <- read.csv("UKB.AMC.current_address.UK_MSOA_regions.csv", header=T, stringsAsFactors = FALSE)
LA.BP <- read.csv("UKB.AMC.birth_place.UK_local_authorities.csv", header=T, stringsAsFactors = FALSE)
MSOA.BP <- read.csv("UKB.AMC.birth_place.UK_MSOA_regions.csv", header=T, stringsAsFactors = FALSE)

colnames(LA.CA)
colnames(MSOA.CA)
colnames(LA.BP)
colnames(MSOA.BP)

LA.CA <- LA.CA[c("IID","geo_code")]
MSOA.CA <- MSOA.CA[c("IID","geo_code")]
LA.BP <- LA.BP[c("IID","geo_code")]
MSOA.BP <- MSOA.BP[c("IID","geo_code")]

colnames(LA.CA) <- c("IID","LA.CA.geo_code")
colnames(MSOA.CA) <- c("IID","MSOA.CA.geo_code")
colnames(LA.BP) <- c("IID","LA.BP.geo_code")
colnames(MSOA.BP) <- c("IID","MSOA.BP.geo_code")

UKB <- merge(UKB.PRS.res, LA.CA, by="IID", all.x=T)
UKB <- merge(UKB, MSOA.CA, by="IID", all.x=T)
UKB <- merge(UKB, LA.BP, by="IID", all.x=T)
UKB <- merge(UKB, MSOA.BP, by="IID", all.x=T)

covar=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_covar_regenie.txt",header=T,stringsAsFactors=FALSE)
covar=as.data.frame(covar)
covar2=covar[match(train_ids,covar$IID),]
write.table(covar2,"/pl/active/KellerLab/Yongkang/PGS_web/UKB_covar_regenie_train.txt",row.names=F,col.names=T,quote=F,sep="\t")

test_ids=fread("/projects/yoki5348/test_samp_ID_all_spo.txt",header=F,stringsAsFactors=FALSE)
test_ids=as.data.frame(test_ids)[,1]
pheno2=pheno[match(test_ids,pheno$IID),]
write.table(pheno2,"/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_test.txt",row.names=F,col.names=T,quote=F,sep="\t")

covar=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_covar_regenie.txt",header=T,stringsAsFactors=FALSE)
covar=as.data.frame(covar)
covar2=covar[match(test_ids,covar$IID),]
write.table(covar2,"/pl/active/KellerLab/Yongkang/PGS_web/UKB_covar_regenie_test.txt",row.names=F,col.names=T,quote=F,sep="\t")
