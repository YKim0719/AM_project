##Systolic blood pressure: 411:418
#Diastolic blood pressure: 419:426
#Birth weight known: 453:455
#Average total household income before tax: 755:758
#Sleeplessness: 910:913
#Snoring: 914:917
#Current tobacco smoking 926:929
#Alcohol intake frequency: 1046:1049
#Mood swing: 1143:1146
#Miserableness 1147:1150
#Irritability 1151:1154
#Sensitivity (hurt feeling): 1155:1158
#Fed-up feelings : 1159-1162
#Nervous feeling: 1163:1166
#Worrier/ anxious feelings: 1167:1170
#Tense (Highly strung): 1171-1174
#Number of cigarettes previously 1445-1448
#FEV1 1537:1548
#Heel bone mineral density 1627
#Number of cigarettes currently smoked daily 1688:1691
#Adopted father still alive: 1839:1842
#Adopted mother still alive: 1843:1846
#Systolic blood pressure (auto): 1887:1894
#Diastolic blood pressure(auto) : 1879:1886
#Bone mineral density (left): 1935:1938
#Bone mineral density (right): 1959:1962
#Fluid intelligence score 9703:9706
#Forced expiratory volume in 1-second (FEV1), Best measure: 10875
#Forced vital capacity (FVC), Best measure: 10876
#Body fat percentage : 13737-13740
#Whole body fat mass: 13741:13744
#Glucose 16916:16917
#HBA1c 16930:16931
#HDL 16942:16943
#LDL 16970:16971
#TG 17096:17097

library(data.table)
library(dplyr)
library(reshape2)

# Import large data file:
load("/pl/active/IBG/UKBiobank/PHENO/Keller/old_baskets/Keller-42049/ukb42049.RData")

FTE=bd[,787:789]#Age completed full time education
BMI=bd[,c(11777:11780)]
#BMI[which(is.na(BMI[,1])),1]=BMI[which(is.na(BMI[,1])),2]
#BMI[which(is.na(BMI[,1])),1]=BMI[which(is.na(BMI[,1])),3]
#BMI[which(is.na(BMI[,1])),1]=BMI[which(is.na(BMI[,1])),4]

PCs=bd[,c(12025:12064)]

total_year=c()
for(check in 1:nrow(FTE)){
  total_year=c(total_year,max(as.numeric(FTE[check,]),na.rm=TRUE))
}

total_year[total_year-5<7]=NA
eduyears=rep(NA,nrow(bd))
for(check in 6981:7004){
  aa=rep(NA,nrow(bd))
  aa[which(bd[,check]=="Prefer not to answer")]=NA
  aa[which(bd[,check]=="College or University degree")]=20
  aa[which(bd[,check]=="A levels/AS levels or equivalent")]=13
  aa[which(bd[,check]=="O levels/GCSEs or equivalent")]=10
  aa[which(bd[,check]=="CSEs or equivalent")]=10
  aa[which(bd[,check]=="NVQ or HND or HNC or equivalent")]=total_year[which(bd[,check]=="NVQ or HND or HNC or equivalent")]-5## Re generate need
  aa[which(bd[,check]=="Other professional qualifications eg: nursing, teaching")]=15
  aa[which(bd[,check]=="None of the above")]=7
  eduyears=cbind(eduyears,aa)
}

eduyears=eduyears[,-1]
aa=c()
for(check in 1:nrow(eduyears)){
  aa=c(aa,max(as.numeric(eduyears[check,]),na.rm=TRUE))
  
}
aa2=aa
aa2[which(aa==-Inf)]=NA
mean(aa2,na.rm=TRUE)

genetic_measure_batch=	bd[,12016] ##genotype measurement batch (UKB field ID: 22000-0.0) 
assessment_center=bd[,c(102:105)]#categorical covariates which were: assessment center (UKB field ID: 54-2.0)

dat=cbind(bd[,c(1,27,28,89,93,102,1131,1465)],aa2,BMI[,1])
dat=cbind(dat,PCs[,1:20],genetic_measure_batch,assessment_center)
colnames(dat)[1:30]=c("IID","Sex","Year_of_birth","Standing_height","Seated_height","UKBiobank_centre","Mother_age","Father_age","EA","BMI",paste("PC",1:20,sep=""))
dat$Year_of_birth2=dat$Year_of_birth^2
dat$UKBiobank_centre=paste("center_",dat$UKBiobank_centre,sep="")
dat$genetic_measure_batch=paste("batch_",dat$genetic_measure_batch,sep="")
write.csv(dat,"/pl/active/KellerLab/Yongkang/UKBhap/EA_modified_pheno_Aysu.csv",row.names=F,quote=F)


###########################################################The other traits ############################################
birth=bd[,453:455]
income=bd[,755]
sleep=bd[,910]
snore=bd[,914]
smoking=bd[,926]
alcohol=bd[,1046]
mood=bd[,1143]
misery=bd[,1147]
irritability=bd[,1151]
sens=bd[,1155]
fed_up=bd[,1159]
nervous=bd[,1163]
worry=bd[,1167]
tense=bd[,1171]
cpd=bd[,1445]
adopt=bd[,1839:1842]
sys=bd[,1887]
dis=bd[,1879]
intel_score=bd[,9703:9706]
fev1=bd[,10875]
fvc=bd[,10876]
body_fat=bd[,13737]
glu=bd[,16916]
hba1c=bd[,16930]
hdl=bd[,16942]
ldl=bd[,16970]
tg=bd[,17096]

dat=bd[,c(1,755,910,914,926,1046,1143,1147,1151,1155,1159,1163,1167,1171,1445,1887,1879,9703,10875,10876,13737,16916,16930,16942,16970,17096)]
colnames(dat)=c("IID","income","sleep","snore","smoking","alcohol","mood","misery","irritability","sensitivity","fed_up","nervous","worry","tense",
                "cpd","adopt","systolic_blood_pressure","diastolic_blood_pressure","intel_flood_score","FEV1","FVC",
                "body_fat","hba1c","hdl","ldl","tg")
dat$income=as.character(dat$income)
aa=names(table(dat$income))
dat$income[which(dat$income==aa[1])]=((18000+31000)/2)/10000
dat$income[which(dat$income==aa[2])]=((31000+52000)/2)/10000
dat$income[which(dat$income==aa[3])]=((52000+100000)/2)/10000
dat$income[which(dat$income==aa[5])]=(100000)/10000
dat$income[which(dat$income==aa[6])]=(18000)/10000
dat$income=as.numeric(dat$income)

dat$sleep=as.character(dat$sleep)
aa=names(table(dat$sleep))
dat$sleep[dat$sleep==aa[1]]=0
dat$sleep[dat$sleep==aa[3]]=1
dat$sleep[dat$sleep==aa[4]]=2
dat$sleep=as.numeric(dat$sleep)

dat$snore=as.character(dat$snore)
aa=names(table(dat$snore))
dat$snore[dat$snore==aa[2]]=1
dat$snore[dat$snore==aa[4]]=2
dat$snore=as.numeric(dat$snore)

dat$smoking=as.character(dat$smoking)
aa=names(table(dat$smoking))
dat$smoking[dat$smoking==aa[1]]=0
dat$smoking[dat$smoking==aa[2]]=1
dat$smoking[dat$smoking==aa[4]]=2
dat$smoking[dat$smoking==aa[5]]=3
dat$smoking=as.numeric(dat$smoking)

dat$alcohol=as.character(dat$alcohol)
aa=names(table(dat$alcohol))
dat$alcohol[dat$alcohol==aa[2]]=0
dat$alcohol[dat$alcohol==aa[6]]=1
dat$alcohol[dat$alcohol==aa[4]]=2
dat$alcohol[dat$alcohol==aa[3]]=3
dat$alcohol[dat$alcohol==aa[7]]=4
dat$alcohol[dat$alcohol==aa[1]]=5
dat$alcohol=as.numeric(dat$alcohol)

dat$mood=as.character(dat$mood)
aa=names(table(dat$mood))
dat$mood[dat$mood==aa[2]]=1
dat$mood[dat$mood==aa[4]]=2
dat$mood=as.numeric(dat$mood)

dat$worry=as.character(dat$worry)
aa=names(table(dat$worry))
dat$worry[dat$worry==aa[2]]=1
dat$worry[dat$worry==aa[4]]=2
dat$worry=as.numeric(dat$worry)

dat$tense=as.character(dat$tense)
aa=names(table(dat$tense))
dat$tense[dat$tense==aa[2]]=1
dat$tense[dat$tense==aa[4]]=2
dat$tense=as.numeric(dat$tense)

dat$misery=as.character(dat$misery)
aa=names(table(dat$misery))
dat$misery[dat$misery==aa[2]]=1
dat$misery[dat$misery==aa[4]]=2
dat$misery=as.numeric(dat$misery)

dat$nervous=as.character(dat$nervous)
aa=names(table(dat$nervous))
dat$nervous[dat$nervous==aa[2]]=1
dat$nervous[dat$nervous==aa[4]]=2
dat$nervous=as.numeric(dat$nervous)

dat$fed_up=as.character(dat$fed_up)
aa=names(table(dat$fed_up))
dat$fed_up[dat$fed_up==aa[2]]=1
dat$fed_up[dat$fed_up==aa[4]]=2
dat$fed_up=as.numeric(dat$fed_up)

dat$sensitivity=as.character(dat$sensitivity)
aa=names(table(dat$sensitivity))
dat$sensitivity[dat$sensitivity==aa[2]]=1
dat$sensitivity[dat$sensitivity==aa[4]]=2
dat$sensitivity=as.numeric(dat$sensitivity)

dat$irritability=as.character(dat$irritability)
aa=names(table(dat$irritability))
dat$irritability[dat$irritability==aa[2]]=1
dat$irritability[dat$irritability==aa[4]]=2
dat$irritability=as.numeric(dat$irritability)
fwrite(dat,"/pl/active/KellerLab/Yongkang/UKBhap/Meng/new_LD/indep_UKB/pheno_all_not_adjusted.csv",row.names=F,quote=F)
pheno=fread("/pl/active/KellerLab/Yongkang/UKBhap/Meng/PGS/210517_NEW/pheno_std_check.csv",header=T,stringsAsFactors=FALSE)
pheno=as.data.frame(pheno)
pheno=pheno[,c(1,2,3,4,9,10,17:116)]
dat2=merge(dat,pheno,by="IID")
write.csv(dat2,"/pl/active/KellerLab/Yongkang/UKBhap/UKBiobank_pheno_all_numeric.csv",row.names=FALSE,quote=FALSE)

colnames(dat2)[1]=c("#IID")
dat3=dat2[,c(1:26,29,30)]
write.table(dat3,"/pl/active/KellerLab/Yongkang/UKBhap/prs_external/UKB_without_used_samp/UKBiobank_phenofile.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)

dat3=dat2[,c(1,27,28,31:131)]
dat3$Year_of_birth2=(dat3$Year_of_birth-mean(dat3$Year_of_birth))^2
write.table(dat3,"/pl/active/KellerLab/Yongkang/UKBhap/prs_external/UKB_without_used_samp/UKBiobank_covfile.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)



#################MDD##########################################

mdd=bd[,17391:17456]

aa=c()
for(idx in 1:ncol(mdd)){
  aa=c(aa,which(mdd[,idx]=="F32"|mdd[,idx]=="F33"|mdd[,idx]=="F34"|mdd[,idx]=="F38"|mdd[,idx]=="F39"))
}
mdd2=bd[,17485:17668]
for(idx in 1:ncol(mdd2)){
  aa=c(aa,which(mdd2[,idx]=="F32"|mdd2[,idx]=="F33"|mdd2[,idx]=="F34"|mdd2[,idx]=="F38"|mdd2[,idx]=="F39"))
}

mdd3=bd[,1211:1218]
for(idx in 1:ncol(mdd3)){
  aa=c(aa,which(mdd3[,idx]=="Yes"))
}

bb=which(mdd3[,1]=="No")
bb2=which(mdd3[,5]=="No")
bb3=which(mdd3[,5]=="No"&mdd3[,1]=="No")
bdd=rep(NA,nrow(mdd))
bdd[bb3]=1
bdd[aa]=2


mdd3=bd[,1211:1218]
aa=c()
for(idx in 1:ncol(mdd3)){
  aa=c(aa,which(mdd3[,idx]=="Yes"))
}
aa=unique(aa)
mdd1=aa[which(aa%in%which(bd[,4099]=="Yes" &bd[,4103]>2))]
mdd2=aa[which(aa%in%which(bd[,4111]=="Yes" &bd[,5417]>2))]
mdds=unique(c(mdd1,mdd2))

prob_mdd=rep(NA,nrow(bd))
mdd=bd[,17391:17456]
aa=c()
for(idx in 1:ncol(mdd)){
  aa=c(aa,which(mdd[,idx]=="F32"|mdd[,idx]=="F33"|mdd[,idx]=="F34"|mdd[,idx]=="F38"|mdd[,idx]=="F39"))
}
mdd2=bd[,17485:17668]
for(idx in 1:ncol(mdd2)){
  aa=c(aa,which(mdd2[,idx]=="F32"|mdd2[,idx]=="F33"|mdd2[,idx]=="F34"|mdd2[,idx]=="F38"|mdd2[,idx]=="F39"))
}
aa=c(aa,mdds)
prob_mdd[aa]=2
##Control
mdd3=bd[,1211:1218]
aa=c()
for(idx in 1:ncol(mdd3)){
  aa=c(aa,which(mdd3[,idx]=="No"))
}
aa=c(aa,which(bd[,4099]=="No"))
aa=c(aa,which(bd[,4111]=="No"))
prob_mdd[aa]=1

icd_mdd=rep(NA,nrow(bd))
mdd=bd[,17391:17456]
aa=c()
for(idx in 1:ncol(mdd)){
  aa=c(aa,which(mdd[,idx]=="F32"|mdd[,idx]=="F33"|mdd[,idx]=="F34"|mdd[,idx]=="F38"|mdd[,idx]=="F39"))
}
mdd2=bd[,17485:17668]
for(idx in 1:ncol(mdd2)){
  aa=c(aa,which(mdd2[,idx]=="F32"|mdd2[,idx]=="F33"|mdd2[,idx]=="F34"|mdd2[,idx]=="F38"|mdd2[,idx]=="F39"))
}
icd_mdd[aa]=2

##Control
mdd3=bd[,1211:1218]
aa=c()
for(idx in 1:ncol(mdd3)){
  aa=c(aa,which(mdd3[,idx]=="No"))
}
aa=c(aa,which(bd[,4099]=="No"))
aa=c(aa,which(bd[,4111]=="No"))
icd_mdd[aa]=1

exclude_icd=c()
mdd=bd[,17391:17456]
for(idx in 1:ncol(mdd)){
  exclude_icd=c(exclude_icd,which(mdd[,idx]=="F30"|mdd[,idx]=="F31"|mdd[,idx]=="F44.8"|mdd[,idx]=="F21"|mdd[,idx]=="F28"|mdd[,idx]=="F29"|mdd[,idx]=="1291"|mdd[,idx]=="1289"))
}
mdd2=bd[,17485:17668]
for(idx in 1:ncol(mdd2)){
  exclude_icd=c(exclude_icd,which(mdd2[,idx]=="F30"|mdd2[,idx]=="F31"|mdd2[,idx]=="F44.8"|mdd2[,idx]=="F21"|mdd2[,idx]=="F28"|mdd2[,idx]=="F29"|mdd2[,idx]=="1291"|mdd2[,idx]=="1289"))
}

treat=bd[,8514:8672]

exclude_tab=c()
for(idx in 1:ncol(treat)){
  exclude_tab=c(exclude_tab,which(treat[,idx]%in%c(1141202024, 1141153490 , 1141195974 , 1140867078 ,1140867494 , 1141171566 , 2038459704 ,1140872064 ,1140879658 ,1140867342 ,1140867420 , 1140882320 , 1140872216 , 1140910358 , 1141200458 , 1141172838 , 1140867306 ,1140867180 ,1140872200,1140867210 ,1140867398 ,1140882098 ,1140867184 ,1140867168 ,1140863416 , 1140909802,1140867498 ,1140867490 ,1140910976 , 1140867118 , 1140867456 , 1140928916 ,1140872268 , 1140867134 ,1140867208 ,1140867218 , 1140867572, 1140879674,1140909804, 1140867504 , 1140868170 , 1140879746,1141152848 ,1141177762, 1140867444 , 1140867092 , 1141152860 , 1140872198 , 1140867244 ,1140868172,1140867304,
                                                   1140872072, 1140879750,1140868120,1140872214, 1141201792,1140882100, 1141167976)))
}

exclude_tab_cont=c()
cont=which(prob_mdd==1)

for(idx in 1:ncol(treat)){
  exclude_tab_cont=c(exclude_tab_cont,cont[cont%in%which(treat[,idx]%in%c(1140867820, 1140867948, 1140879616, 1140867938, 1140867690,1141190158,1141151946,1140921600,1140879620,1141201834,1140867152,
                                                                          1140909806, 1140879628, 1140867640, 1141200564, 1141151982, 1140916288, 1141180212, 1140867860, 1140867952, 1140879540,
                                                                          1140867150, 1140909800, 1140867940, 1140879544, 1140879630, 1140867856, 1140867726, 1140867884, 1140867922, 1140910820,
                                                                          1140879556, 1141152732, 1140867920, 1140882244, 1140867852, 1140867818, 1141174756, 1140867916, 1140867888, 1140867850,
                                                                          1140867624, 1140867876, 1141151978, 1140882236, 1140867878, 1201, 1140882312, 1140867758, 1140867712, 1140867914, 1140867944,
                                                                          1140879634, 1140867756, 1140867934, 1140867960, 1140916282, 1141200570, 1141152736))])
}

mdd=bd[,17391:17456]
aa=c()
for(idx in 1:ncol(mdd)){
  aa=c(aa,which(mdd[,idx]=="F32"|mdd[,idx]=="F33"|mdd[,idx]=="F34"|mdd[,idx]=="F38"|mdd[,idx]=="F39"|mdd[,idx]=="1286"))
}
mdd2=bd[,17485:17668]
for(idx in 1:ncol(mdd2)){
  aa=c(aa,which(mdd2[,idx]=="F32"|mdd2[,idx]=="F33"|mdd2[,idx]=="F34"|mdd2[,idx]=="F38"|mdd2[,idx]=="F39"|mdd2[,idx]=="1286"))
}

exclude_tab_cont=c(exclude_tab_cont,cont[cont%in%aa])

prob_mdd[exclude_tab]=NA
bdd[exclude_tab]=NA
icd_mdd[exclude_tab]=NA
prob_mdd[exclude_tab_cont]=NA
bdd[exclude_tab_cont]=NA
icd_mdd[exclude_tab_cont]=NA
prob_mdd[exclude_icd]=NA
bdd[exclude_icd]=NA
icd_mdd[exclude_icd]=NA

dat=bd[,c(1,27,28,89,102)]
colnames(dat)=c("IID","Sex","Year_of_birth","Standing_height","UKBiobank_centre")
dat=cbind(dat,bdd,prob_mdd,icd_mdd)

dat2=fread("/pl/active/KellerLab/Yongkang/UKBhap/prs_external/UKB_without_used_samp/UKBiobank_phenofile_adjusted.txt",header=T,stringsAsFactors=FALSE)
dat2=as.data.frame(dat2)
dat_european=dat[match(dat2$"#IID",dat$IID),]
dat_european=dat_european[,-c(2,3,5)]
dat_european=dat_european[,-2]
colnames(dat_european)[1]="#IID"
write.table(dat_european,"/pl/active/KellerLab/Yongkang/PGS_web/UKB_phenofile_mdd.txt",row.names=F,col.names=T,quote=F,sep="\t")





aa=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_phenofile_mdd.txt",header=T,stringsAsFactors=FALSE)
aa=as.data.frame(aa)
cor(pheno2$pgs_major_depression_non_UKB,aa[,2],use="complete.obs")
cor(pheno2$pgs_major_depression_non_UKB,aa[,3],use="complete.obs")
cor(pheno2$pgs_major_depression_non_UKB,aa[,4],use="complete.obs")


aa=fread("/pl/active/KellerLab/Yongkang/PGS_web/pheno_richard_mdd.txt",header=T,stringsAsFactors=FALSE)
aa=as.data.frame(aa)
aa=aa[match(pheno2$IID,aa$id),]
cor(pheno2$pgs_major_depression_non_UKB,aa[,2],use="complete.obs")
cor(pheno2$pgs_major_depression_non_UKB,aa[,3],use="complete.obs")
cor(pheno2$pgs_major_depression_non_UKB,aa[,4],use="complete.obs")

aa=dir("/pl/active/KellerLab/Yongkang/PGS_web/ldsc/intelligence/whole_variant/individual_hpgs2")
bb=paste("chr",c(rep(1:12,each=50),rep(13:22,each=10)),".",c(rep(1:50,12),rep(1:10,10)),"_ver2.csv",sep="")
bb[-which(bb%in%aa)]


#####life time MDD ####
#To receive a diagnosis, participants had to meet all of the following criteria:
#1.Anhedonia or depressed mood. Respondents needed to have responded affirmatively to either of the following questions
##Anhedonia (20441) "Have you ever had a time in your life lasting two weeks or more when you lost interest in most things like hobbies, work, or activities that usually give you pleasure?"
anhedonia=bd[,11477]
##Mood (20446) "Have you ever had a time in your life when you felt sad, blue, or depressed for two weeks or more in a row?"
mood=bd[,11480]

lifetime1=rep(NA,nrow(bd))
lifetime1[which(anhedonia=="Yes")]=2
lifetime1[which(anhedonia=="No")]=1
lifetime2=rep(NA,nrow(bd))
lifetime2[which(mood=="Yes")]=2
lifetime2[which(mood=="No")]=1

lifetime=rep(NA,nrow(bd))
lifetime[which(lifetime1==2|lifetime2==2)]=2
lifetime[which(lifetime1==1&lifetime2==1)]=1


#2.Symptom count. Respondents needed to endorse 4 or more of the following symptoms (unfortunately, motor agitation/retardation was not assessed) with respect to their worst period of depression
##Anhedonia (20441) "Have you ever had a time in your life lasting two weeks or more when you lost interest in most things like hobbies, work, or activities that usually give you pleasure?"
anhedonia1=bd[,11477]
anhedonia=rep(1,nrow(bd))
anhedonia[anhedonia1=="Yes"]=2
##Mood (20446) "Have you ever had a time in your life when you felt sad, blue, or depressed for two weeks or more in a row?"
mood1=bd[,11480]
mood=rep(1,nrow(bd))
mood[mood1=="Yes"]=2
##Sleep (20533, 20534, 20533) “Trouble falling asleep” or “sleeping too much” or “waking too early”
sleep1=bd[,11560]
sleep2=bd[,11561]
sleep3=bd[,11562]
sleep=rep(1,nrow(bd))
sleep[which(sleep1=="Yes"|sleep2=="Yes"|sleep3=="Yes")]=2
##Fatigue (20449) "Did you feel more tired out or low on energy than is usual for you?"
fatigue1=bd[,11483]
fatigue=rep(1,nrow(bd))
fatigue[fatigue1=="Yes"]=2
##Appetite/weight (20536) "Did you gain or lose weight without trying, or did you stay about the same weight?
appetite1=bd[,11563]
appetite=rep(1,nrow(bd))
appetite[appetite1=="Gained weight"|appetite1=="Lost weight"|appetite1=="Both gained and lost some weight during the episode"]=2
##Feelings of worthlessness (20450) "People sometimes feel down on themselves, no good, worthless. Did you feel this way?"
worthless1=bd[,11484]
worthless=rep(1,nrow(bd))
worthless[worthless1=="Yes"]=2
##Concentration (20435) "Did you have a lot more trouble concentrating than usual?"
concentration1=bd[,11471]
concentration=rep(1,nrow(bd))
concentration[concentration1=="Yes"]=2
##Ideation (20437) "Did you think a lot about death - either your own, someone else’s, or death in general?"
ideation1=bd[,11473]
ideation=rep(1,nrow(bd))
ideation[ideation1=="Yes"]=2
conti_symp=anhedonia+mood+sleep+fatigue+appetite+worthless+concentration+ideation-8
symptom=as.numeric(anhedonia+mood+sleep+fatigue+appetite+worthless+concentration+ideation-8>=4)+1

#3.Frequency. (20439) With respect to their worst period of depression, respondents had to indicate a 2 or 3 on the following scale assessing “how often [they] felt this way”: {1: “Less often”, 2: “Almost every day”, 3: “Every day”}
frequency1=bd[,11475]
frequency=as.numeric(frequency1=="Almost every day"|frequency1=="Every day")+1

#4.Fraction of day affected. (20436) With respect to their worst period of depression, respondents had to indicate a 2 or greater on the following scale assessing “How much of the day did these feelings usually last?”: {1: “Less than half of the day”, 2: “About half of the day”, 3: “Most of the day”, 4: “All day long”}
fraction1=bd[,11472]
fraction=as.numeric(fraction1=="About half of the day"|fraction1=="Most of the day"| fraction1=="All day long")+1
#5.Impairment. (20440) With respect to their worst period of depression, respondents had to indicate a 2 or greater on the following scale assessing "Think about your roles at the time of this episode, including study / employment, childcare and housework, leisure pursuits. How much did these problems interfere with your life or activities?" : {0: “Not at all”, 1: “A little”, 2: “Somewhat”, 3: “A lot”}
impairment1=bd[,11476]
impairment=as.numeric(impairment1=="Somewhat"|impairment1=="A lot")+1

#exckude
mania=bd[,17391:17456]
aa=c()
for(idx in 1:ncol(mania)){
  aa=c(aa,which(mania[,idx]=="R440"|mania[,idx]=="R441"|mania[,idx]=="R442"|mania[,idx]=="R443"))
  aa=c(aa,which(mania[,idx]=="F220"|mania[,idx]=="F22"|mania[,idx]=="F24"|mania[,idx]=="F223"|mania[,idx]=="F228"|mania[,idx]=="F229"))
}

mania2=bd[,17485:17668]
for(idx in 1:ncol(mania)){
  aa=c(aa,which(mania2[,idx]=="R440"|mania2[,idx]=="R441"|mania2[,idx]=="R442"|mania2[,idx]=="R443"))
  aa=c(aa,which(mania2[,idx]=="F220"|mania2[,idx]=="F22"|mania2[,idx]=="F24"|mania2[,idx]=="F223"|mania2[,idx]=="F228"|mania2[,idx]=="F229"))
}

life_mdd=as.numeric(lifetime==2&symptom==2&frequency==2&fraction==2&impairment==2)

##### Current MDD severity #####
# Respondents ranked the severity of each symptom on a 0-4 scale rather than providing a binary endorsement
# All symptoms were assessed regardless of whether or not anhedonia or depressed mood was diagnosed.
#The following items were assessed as “Over the last 2 weeks, how often have you been bothered by any of the following problems?”
## 1. Anhedonia (20510) "Little interest or pleasure in doing things"
anhedonia1=bd[,11537]
anhedonia=rep(NA,nrow(bd))
anhedonia[anhedonia1=="Not at all"]=0
anhedonia[anhedonia1=="Several days"]=1
anhedonia[anhedonia1=="More than half the days"]=2
anhedonia[anhedonia1=="Nearly every day"]=3

##2. Mood (20514) "Feeling down, depressed, or hopeless"
mood1=bd[,11541]
mood=rep(NA,nrow(bd))
mood[mood1=="Not at all"]=0
mood[mood1=="Several days"]=1
mood[mood1=="More than half the days"]=2
mood[mood1=="Nearly every day"]=3

#3. Sleep (20517) “Trouble falling or staying asleep, or sleeping too much”
sleep1=bd[,11544]
sleep=rep(NA,nrow(bd))
sleep[sleep1=="Not at all"]=0
sleep[sleep1=="Several days"]=1
sleep[sleep1=="More than half the days"]=2
sleep[sleep1=="Nearly every day"]=3

#4. Fatigue (20519) "Did you feel more tired out or low on energy than is usual for you?"
fatigue1=bd[,11546]
fatigue=rep(NA,nrow(bd))
fatigue[fatigue1=="Not at all"]=0
fatigue[fatigue1=="Several days"]=1
fatigue[fatigue1=="More than half the days"]=2
fatigue[fatigue1=="Nearly every day"]=3

#5. Appetite/weight (20511) "Poor appetite or overeating"
appetite1=bd[,11538]
appetite=rep(NA,nrow(bd))
appetite[appetite1=="Not at all"]=0
appetite[appetite1=="Several days"]=1
appetite[appetite1=="More than half the days"]=2
appetite[appetite1=="Nearly every day"]=3

#6. Feelings of worthlessness (20507) "Feeling bad about yourself or that you are a failure or have let yourself or your family down"
worthless1=bd[,11534]
worthless=rep(NA,nrow(bd))
worthless[worthless1=="Not at all"]=0
worthless[worthless1=="Several days"]=1
worthless[worthless1=="More than half the days"]=2
worthless[worthless1=="Nearly every day"]=3

#7. Concentration (20508) "Trouble concentrating on things, such as reading the newspaper or watching television"
concentration1=bd[,11535]
concentration=rep(NA,nrow(bd))
concentration[concentration1=="Not at all"]=0
concentration[concentration1=="Several days"]=1
concentration[concentration1=="More than half the days"]=2
concentration[concentration1=="Nearly every day"]=3

#8. Ideation (20513) "Thoughts that you would be better off dead or of hurting yourself in some way"
ideation1=bd[,11540]
ideation=rep(NA,nrow(bd))
ideation[ideation1=="Not at all"]=0
ideation[ideation1=="Several days"]=1
ideation[ideation1=="More than half the days"]=2
ideation[ideation1=="Nearly every day"]=3

mdd_severe=anhedonia+mood+sleep+fatigue+worthless+appetite+concentration+ideation
mdd_severe[aa]=NA


#3.1.4 Lifetime episode count
#Ordinal measure of lifetime number of depressive episodes among UKBB online mental health follow-up respondents
#who endorsed anhedonia and/or depressed mood (20442).
#Individuals who endorsed a two week period of either anhedonia or depressed mood were asked "How many periods
#did you have in your life lasting two or more weeks where you felt like this?". Respondents supplied either an integer
#between 1 and 999 or responded “Too many to count / One episode ran into the next”, rendering counts greater
#than one difficult to compare. We thus assigned scores as follows:
#  0 individuals who endorsed neither anhedonia or depressed mood
#1 individuals who indicated a single depressive episode
#2 individuals who indicated ≥ 2 depressive episodes

pheno=data.frame(id=bd[,1],life_mdd=life_mdd,conti_symp=conti_symp,mdd_severe=mdd_severe)





write.table(pheno,"/pl/active/KellerLab/Yongkang/PGS_web/pheno_richard_mdd.txt",row.names=F,col.names=T,quote=F)




########################################
sys=bd[,1887]
dis=bd[,1879]
intel_score=bd[,9703:9706]
fev1=bd[,10875]
fvc=bd[,10876]
body_fat=bd[,13737]
glu=bd[,16916]
hba1c=bd[,16930]
hdl=bd[,16942]
ldl=bd[,16970]
tg=bd[,17096]

sleepdur=bd[,894]
gettingup=bd[,898]
morningness=bd[,902]
napping=bd[,906]
insomnia=bd[,910]
snoring=bd[,914]
dozing=bd[,918]

#Field ID	Description
#1160	Sleep duration
#1170	Getting up in morning
#1180	Morning/evening person (chronotype)
#1190	Nap during day
#1200	Sleeplessness / insomnia
#1210	Snoring
#1220	Daytime dozing / sleeping (narcolepsy)


smoke_init_prev=bd[,1437] 
smoke_init_prev[smoke_init_prev==-3]=NA
smoke_init_prev[smoke_init_prev==-1]=NA
smoke_init_current=bd[,1680]
smoke_init_current[smoke_init_current==-3]=NA
smoke_init_current[smoke_init_current==-1]=NA
AgeSmk=smoke_init_prev
AgeSmk[!is.na(smoke_init_current)]=smoke_init_current[!is.na(smoke_init_current)]
cpd=bd[,1445]
cpd[cpd==-10]=NA
cpd[cpd==-1]=NA
cpd2=cpd

cpd2[cpd>=36]=5
cpd2[cpd>=26&cpd<=35]=4
cpd2[cpd>=16&cpd<=25]=3
cpd2[cpd>=6&cpd<=15]=2
cpd2[cpd>=1&cpd<=5]=1
CigDay=cpd2
cpd3=cpd2
cpd3[is.na(cpd2)]=0


ever_smoked=bd[,10908]
current_smoke=bd[,922]
previous_smoke=bd[,926]

SmkCes=rep(NA,nrow(bd))
SmkCes[previous_smoke=="Smoked on most or all days"|
         previous_smoke=="Smoked occasionally"|
         previous_smoke=="Just tried once or twice"]=0
SmkCes[which(current_smoke=="Yes, on most or all days"|current_smoke=="Only occasionally")]=1

SmkInit=rep(NA,nrow(bd))
SmkInit[which(current_smoke=="Yes, on most or all days"|current_smoke=="Only occasionally"|previous_smoke=="Smoked on most or all days"|
                previous_smoke=="Smoked occasionally")]=1
SmkInit[which(previous_smoke=="I have never smoked")]=0

alcohol=bd[,11452]
alcohol2=rep(NA,nrow(bd))
alcohol2[which(alcohol=="Monthly or less")]=0.25
alcohol2[which(alcohol=="2 to 4 times a month")]=0.75
alcohol2[which(alcohol=="2 to 3 times a week")]=2.5
alcohol2[which(alcohol=="4 or more times a week")]=5

alcohol3=rep(NA,nrow(bd))
alcohol3[which(alcohol=="Never")]=0
alcohol3[which(alcohol=="Monthly or less")]=0.25
alcohol3[which(alcohol=="2 to 4 times a month")]=0.75
alcohol3[which(alcohol=="2 to 3 times a week")]=2.5
alcohol3[which(alcohol=="4 or more times a week")]=5

freq_alcohol=bd[,11441]
freq_alcohol2=rep(NA,nrow(bd))
freq_alcohol2[which(freq_alcohol=="1 or 2")]=1.5
freq_alcohol2[which(freq_alcohol=="3 or 4")]=3.5
freq_alcohol2[which(freq_alcohol=="5 or 6")]=5.5
freq_alcohol2[which(freq_alcohol=="7, 8 or 9")]=8
freq_alcohol2[which(freq_alcohol=="10 or more")]=11

freq_alcohol3=rep(0,nrow(bd))
freq_alcohol3[which(freq_alcohol=="1 or 2")]=1.5
freq_alcohol3[which(freq_alcohol=="3 or 4")]=3.5
freq_alcohol3[which(freq_alcohol=="5 or 6")]=5.5
freq_alcohol3[which(freq_alcohol=="7, 8 or 9")]=8
freq_alcohol3[which(freq_alcohol=="10 or more")]=11

DrnkWk=alcohol2*freq_alcohol2
DrnkWk2=alcohol3*freq_alcohol3

dat=data.frame(IID=bd[,1],AgeSmk=AgeSmk,CigDay=CigDay,CigDay_include_never_smoke=cpd3,SmkInit=SmkInit,DrnkWk=DrnkWk,SmkCes=SmkCes,DrnkWk2=DrnkWk2)
write.csv(dat,"/pl/active/KellerLab/Yongkang/PGS_web/UKB_phenofiles_smoke_alcohol.csv",row.names=F,quote=F)


