import os

pair2Info = {}
f = open(snakemake.input[0])
header = f.readline()
for line in f:
    lines = line.strip().split("\t")
    pair = "_".join(lines[:2])
    pValue = float(lines[9])
    if pair not in pair2Info:
        pair2Info[pair] = {"pValue":pValue,"info":lines}
    elif pair in pair2Info:
        if pValue < pair2Info[pair]["pValue"]:
            print(f"{pair}:\t{pair2Info[pair]['pValue']}")
            pair2Info[pair] = {"pValue":pValue,"info":lines}
f.close()


g = open(snakemake.output[0],"w")
g.write(header)
for pair in pair2Info:
    g.write("\t".join(pair2Info[pair]["info"]) + "\n")
g.close()