

library(data.table)

args <- commandArgs(trailingOnly =TRUE)
chr=as.numeric(args[1])
pheno_idx=as.numeric(args[2])

pheno=fread("/pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie.txt",header=T,stringsAsFactors=FALSE)
pheno=as.data.frame(pheno)
dir.create("/pl/active/KellerLab/Yongkang/PGS_web/regenie_train_test")
  phenoid=colnames(pheno)[pheno_idx+2]
  dir.create(paste("/pl/active/KellerLab/Yongkang/PGS_web/regenie_train_test/pheno",phenoid,sep=""))
    system(paste("/pl/active/KellerLab/opt/bin/regenie --bsize 100 --bed /pl/active/KellerLab/Yongkang/hapmap3/all_eur_bed/chr",chr,"_common --strict --phenoCol ",phenoid," --threads 5 --phenoFile /pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_train.txt --covarFile /pl/active/KellerLab/Yongkang/PGS_web/UKB_covar_regenie_train.txt  --pred /pl/active/KellerLab/Yongkang/PGS_web/regenie_train_test/pheno",phenoid,"/chr",chr,".list --step 1 --out /pl/active/KellerLab/Yongkang/PGS_web/regenie_train_test/pheno",phenoid,"/chr",chr,sep="")) 
    system(paste("/pl/active/KellerLab/opt/bin/regenie --bsize 200 --bed /pl/active/KellerLab/Yongkang/hapmap3/all_eur_bed/chr",chr,"_common --strict --phenoCol ",phenoid," --threads 3 --sample /projects/yoki5348/train_samp_ID_w.o.spo.txt --phenoFile /pl/active/KellerLab/Yongkang/PGS_web/UKB_pheno_regenie_train.txt --covarFile /pl/active/KellerLab/Yongkang/PGS_web/UKB_covar_regenie_train.txt  --step 2  --firth --approx --pThresh 0.01 --pred /pl/active/KellerLab/Yongkang/PGS_web/regenie_train_test/pheno",phenoid,"/chr",chr,"_pred.list --out /pl/active/KellerLab/Yongkang/PGS_web/regenie_train_test/pheno",phenoid,"/chr",chr,sep="")) 

