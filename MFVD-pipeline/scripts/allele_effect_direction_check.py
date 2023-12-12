import os
import gzip

f = gzip.open(snakemake.input[0],"rt")
g = open(snakemake.output[0],"w")
g.write("\t".join(["CHROM","POS","REF","ALT","heterozygosis_sample_count","REF_count","ALT_count","DP_count"]) + "\n")

for line in f:
    if not line.startswith("#"):
        lines = line.strip().split("\t")
        chromosome,pos = lines[:2]
        REF,ALT = lines[3:5]

        heterozygosis_count,REF_count,ALT_count,DP_count,not_equal_count,ref_allele_smaller_count = 0,0,0,0
        
        for i in lines[9:]:
            FORMAT_list = i.split(":")
            if FORMAT_list[0] == "0/1":
                heterozygosis_count += 1
                REF_count += int(FORMAT_list[1].split(",")[0])
                ALT_count += int(FORMAT_list[1].split(",")[1])
                DP_count += int(FORMAT_list[2])

                if int(FORMAT_list[1].split(",")[0]) != int(FORMAT_list[1].split(",")[1]):
                    not_equal_count += 1
                    if int(FORMAT_list[1].split(",")[0]) < int(FORMAT_list[1].split(",")[1]):
                        ref_allele_smaller_count += 1

        if ref_allele_smaller_count == 0 or ref_allele_smaller_count == not_equal_count:
            g.write("\t".join([chromosome,pos,REF,ALT,str(heterozygosis_count),str(REF_count),str(ALT_count),str(DP_count)]) + "\n")
f.close()
g.close()
