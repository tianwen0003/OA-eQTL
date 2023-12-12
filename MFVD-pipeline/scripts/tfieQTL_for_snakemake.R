# conda activate cieQTL

options(stringsAsFactors = FALSE);
options(warn = -1)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(stringr)


#----set input file dir and input file----#
eVariant2TF = snakemake@input[["eVariant2TF"]]
exp = snakemake@input[["exp"]]
cov = snakemake@input[["cov"]]
samplelist = snakemake@input[["samplelist"]]

batchNo <- unlist(strsplit(eVariant2TF,split="_"))[length(unlist(strsplit(eVariant2TF,split="_")))]

#----set output file----#
out = snakemake@output[[1]]

#----process expression data(phenotyp)----#
exp <- read.table(exp, header=T, row.names=4, check.names=F,comment.char="")
exp <- exp[,4:ncol(exp)]

#----process covariant data----#
cov = read.table(cov, header=T, row.names=1, check.names=F)
cov = data.frame(t(cov))
cov = cov[colnames(exp),]

#----calculate tfieQTL----#
eVariant2TF = read.table(eVariant2TF, header=T)

#----get beta, std.Error_of_beta, pvalue_of_beat, pvalue_of_the_diff_of_the_two_model----#
for (i in 1:dim(eVariant2TF)[1]){
    chr = eVariant2TF[i,1] # chromosome
    snp = eVariant2TF[i,4] # rs number
    tf = eVariant2TF[i,16] # tf gene id
    target = eVariant2TF[i,5] # eQTL gene id

    tf_exp = t(exp[tf,])
    target_exp = t(exp[target,])

    cov = cbind(cov,tf_exp)

    formula1 = paste0("target_exp2~tpmdos2+(1|sex)+",paste(colnames(cov)[-which(colnames(cov)=="sex")], collapse = "+"))
    formula2 = paste0(str_c("target_exp2~",tf,":tpmdos2+tpmdos2+(1|sex)+"),paste(colnames(cov)[-which(colnames(cov)=="sex")], collapse = "+"))
    
    write.table(snp,str_c(tf,"_",target,"_",batchNo,".snps.list"),quote=F,col.names = F,row.names = F)
    plink_out = paste0(tf, "_", target, "_", batchNo)
    geno = paste0(plink_out,".raw")
    system(paste("/home/TW/anaconda3/bin/plink --bfile /ext3/TW/project/OA_eQTL/eQTL-cis/01_eqtlInputData/OA.TOP.imputed.merged.filtered.DC", "--keep", samplelist, "--out", plink_out, "--recode A", "--extract", str_c(tf,"_",target,"_",batchNo,".snps.list"), collape = " "))
    geno = read.table(geno, header=T, row.names=1,check.names=F)
    tpmdos = geno[,6]
    if (sum(is.na(tpmdos))!=0){
        target_exp2 = as.matrix(target_exp[-which(is.na(tpmdos)),])
        colnames(target_exp2) = target
        cov2 = cov[-which(is.na(tpmdos)),]
        tpmdos2 = tpmdos[-which(is.na(tpmdos))]
    }else{
        target_exp2 = target_exp
        cov2 = cov
        tpmdos2 = tpmdos
    }

    snp = strsplit(colnames(geno)[6], split="_")[[1]]
    snp = paste(snp[-length(snp)], collapse="_")

    print(tf)
    print(target)
    print(snp)

    tpmmodel2<-lmer(formula2, data = cov2, REML = FALSE)
    tpmmodel1<-lmer(formula1, data = cov2, REML = FALSE)
    
    anova(tpmmodel2,tpmmodel1)->comp1
    comp1$AIC[2]-comp1$AIC[1]->sdiff
    
    compare<-PBmodcomp(tpmmodel2,tpmmodel1,nsim=1000)
    summary(compare)$test$p.value[1]->bpval
    beta = summary(tpmmodel2)$coefficients[str_c(tf,":tpmdos2"),]
    pvalues<-data.frame(chiSquarePvalue=bpval)
    diff<-data.frame(AICDiff=sdiff)
    pair<-data.frame(TF=as.character(tf),target_gene=as.character(target),SNP=as.character(snp))
    inter <- t(beta)
    
    cbind(pair,diff,pvalues, inter)->res
    write.table(res, out, sep="\t", quote=F, col.names=F, row.names=F, append=T)
    system(str_c("rm ",tf,"_",target,"_",batchNo,"*"))
}

#----add header for result file----#
system(paste0("sed -i '1i TF\ttaget_gene\tSNP\tAICDiff\tchiSquarePvalue\tEstimate\tSE\tdf\tt_value\tinteractP' ",out))
print("Done!")
