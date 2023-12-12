rm(list = ls())
library(stringr)

vcf_summary <- snakemake@input[["vcf_summary"]]

a <- read.table(file = vcf_summary, sep="\t", header=T)

all_test <- vector()
for (i in 1:nrow(a)){
  if (a[i,8] >= 20 & a[i,1] != "chrM" & a[i,1] != "chrX" & a[i,1] != "chrY"){
    binom_Test <- binom.test(a[i,7],a[i,8],0.5)
    if (binom_Test$p.value <= 1){
      result <- c(data.frame(t(a[i,]))[,1],binom_Test$p.value)
      all_test <- rbind(all_test,result)
    }
  }
}

colnames(all_test) <- c("CHROM","POS","REF","ALT","heterozygosis_sample_count","REF_count","ALT_count","DP_count","P_value")

all_test <- as.data.frame(all_test)
all_test$heterozygosis_sample_count <- as.numeric(all_test$heterozygosis_sample_count)
all_test$REF_count <- as.numeric(all_test$REF_count)
all_test$ALT_count <- as.numeric(all_test$ALT_count)
all_test$DP_count <- as.numeric(all_test$DP_count)
all_test$P_value <- as.numeric(all_test$P_value)
all_test$REF_count <- as.numeric(all_test$REF_count)

all_test$FDR <- p.adjust(all_test$P_value,method = "fdr")
all_test <- all_test[order(all_test$P_value),]

write.table(all_test,file = snakemake@output[["all_"]],quote = F,row.names = F,col.names = T,sep = "\t")
write.table(significant_pvalue,file = snakemake@output[["p"]],quote = F,row.names = F,col.names = T,sep = "\t")
write.table(significant_fdr_0.1,file = snakemake@output[["fdr_1"]],quote = F,row.names = F,col.names = T,sep = "\t")
write.table(significant_fdr_0.05,file = snakemake@output[["fdr_2"]],quote = F,row.names = F,col.names = T,sep = "\t")

rm(list = ls())