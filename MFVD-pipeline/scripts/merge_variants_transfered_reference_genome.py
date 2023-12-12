import os

chrList = ["chr"+str(i) for i in range(1,23)]
chrList += ["chrX","chrY","chrM"]

g = open(snakemake.output[0],"w")
for chr_ in chrList:
    g.write(f">{chr_}\n")
    f = open(f"{snakemake.params}/{chr_}.SNPs_introduced.fa")
    f.readline()
    for line in f:
        g.write(line)
    f.close()
g.close()