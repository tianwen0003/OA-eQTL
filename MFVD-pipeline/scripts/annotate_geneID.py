import os

geneid2Name = {}
f = open(snakemake.input[0])
for line in f:
    lines = line.strip().split("\t")
    geneid2Name[lines[0]] = lines[1]
f.close()

f = open(snakemake.input[1])
g = open(snakemake.output[0],"w")
g.write("TF_id\tTF_name\ttarget_id\ttarget_name\t" + "\t".join(f.readline().strip().split("\t")[2:]) + "\n")
for line in f:
    lines = line.strip().split("\t")
    TF_id,target_id = lines[:2]
    TF_name = geneid2Name[TF_id]
    target_name = geneid2Name[target_id]
    g.write("\t".join([TF_id,TF_name,target_id,target_name]) + "\t" + "\t".join(lines[2:]) + "\n")
f.close()
g.close()
