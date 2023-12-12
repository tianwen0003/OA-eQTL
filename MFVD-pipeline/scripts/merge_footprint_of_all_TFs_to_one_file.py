import os
import re

inpout_dir = snakemake.params[0]
output_dir = snakemake.params[1]

#----load gene name annotation----#
geneName2id = {}
f = open(snakemake.input[2])
for line in f:
    lines = line.strip().split("\t")
    geneName2id[lines[1]] = lines[0]
f.close()

#----merge all footprint result to one file----#
for tfDatabase in ["HOCOMOCO","jaspar"]:
    g = open(f"{output_dir}/{tfDatabase}.all_motif_footprint.bed","w")
    input_dir = f"{inpout_dir}/03_BINDetect_{tfDatabase}_combined"
    
    for file in os.listdir(input_dir):
        if file.endswith("_bound.bed"):
            motif = re.sub("_all_samples_footprints_bound.bed","",file)
            tf_name = motif.split("_")[0]
            f = open(f"{input_dir}/{file}")
            for line in f:
                lines = line.strip().split("\t")
                motif_chr,motif_start,motif_end = lines[:3]
                motif_start = str(int(motif_start)+1)
                motif_strand = lines[5]
                g.write("\t".join([motif_chr,motif_start,motif_end,motif,tf_name,geneName2id[tf_name],motif_strand]) + "\n")
            f.close()
    g.close()