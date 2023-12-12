import os
import gzip

inputFile = snakemake.input[0]
outputFile = snakemake.output[0]

f = gzip.open(inputFile,"rt")
g = gzip.open(outputFile,"wt")

for line in f:
    if line.startswith("##"):
        g.write(line)
    elif line.startswith("#CHROM"):
        lines = line.rstrip("\n").split("\t")
        g.write("\t".join(lines) + "\n")
    else:
        lines = line.rstrip("\n").split("\t")
        formatField = lines[-1].split(":")
        if formatField[0] == "0/1":
            g.write(line)
f.close()
g.close()

tabix_path = snakemake.params[0]
os.system(f"{tabix_path} {outputFile}")