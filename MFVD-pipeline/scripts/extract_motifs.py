import os

#----load expressed genes----#
expressedGenes = set()
f = open(snakemake.input[2])
f.readline()
for line in f:
    lines = line.strip().split("\t")
    if float(lines[1]) >= float(snakemake.params[0]):
        expressedGenes.add(lines[0])
f.close()

#----extract motifs of exprssed TFs----#
for index in [0,1]:
    g = open(snakemake.output[index],"w")
    f = open(snakemake.input[index])
    allLines = f.readlines()
    for i in range(len(allLines)):
        if allLines[i].startswith(">"):
            gene = allLines[i].rstrip("\n").split()[-1]
            if gene in expressedGenes:
                g.write("".join(allLines[i:i+5]))
    f.close()
    g.close()
