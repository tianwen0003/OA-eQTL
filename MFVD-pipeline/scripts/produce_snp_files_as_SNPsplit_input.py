import os
from multiprocessing import Pool

def calibrateRefAllele(line):
    new_lines = []
    lines = line.strip().split("\t")
    chr_= lines[0]
    pos,A1,A2=lines[3:]
    if len(A1) == 1 and len(A2) == 1:
        ref_allele = os.popen(f"samtools faidx {ref_seq} {chr_}:{pos}-{pos} | tail -n 1").read().strip().upper()
        if A1 == ref_allele:
            new_lines = lines[:4] + [A1,A2]
        elif A2 == ref_allele:
            new_lines = lines[:4] + [A2,A1]
    return new_lines

lineList = []
f = open(snakemake.input[0])
for line in f:
    lineList.append(line)
f.close()

threads = int(snakemake.threads)
ref_seq = snakemake.input[1]

with Pool(threads) as pool:
    a = pool.map(calibrateRefAllele,lineList)

chr2snp = {}
for i in a:
    if len(i) > 0:
        chr_,rsID = i[:2]
        chr2snp.setdefault(chr_,[])
        pos,Ref,Alt=i[3:]
        chr2snp[chr_].append([rsID,chr_[3:],pos,"1",Ref+"/"+Alt])

output_dir = snakemake.params

for chr_ in chr2snp:
    print(chr_)
    g = open(f"{output_dir}/{chr_}.txt","w")
    g.write("SNP-ID\tChromosome\tPosition\tStrand\tRef/SNP\n")
    for snp in chr2snp[chr_]:
        g.write("\t".join(snp) + "\n")
    g.close()