import os

output_dir = snakemake.params[0]

f = open(snakemake.input[0])
header = f.readline()
pairs = f.readlines()
f.close()

k,m = divmod(len(pairs),100)
batch = 0

for i in range(100):
    batch += 1
    g = open(f"{output_dir}/variant.TF_motif.overlap_{batch}","w")
    g.write(header)
    for pair in pairs[i*k+min(i, m):(i+1)*k+min(i+1, m)]:
        g.write(pair)
    g.close()