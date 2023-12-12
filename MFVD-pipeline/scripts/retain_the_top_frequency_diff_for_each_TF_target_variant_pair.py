import os

work_dir = snakemake.params[0]
pair2info = {}
for motifdatabase in ["HOCOMOCO","jaspar"]:
    f = open(f"{work_dir}/variant.TF_footprint.overlap.{motifdatabase}.motif_disruption.match_threshold")
    header = f.readline()
    for line in f:
        lines = line.strip().split("\t")
        
        tf_id = lines[15]
        target_id = lines[4]
        variant = lines[3]

        frequecy_diff = float(lines[-1])
        
        pair = "_".join([tf_id,target_id,variant])

        if pair not in pair2info:
            pair2info[pair] = {"info":lines,"frequecy_diff":frequecy_diff}
        elif pair in pair2info:
            if frequecy_diff > pair2info[pair]["frequecy_diff"]:
                pair2info[pair] = {"info":lines,"frequecy_diff":frequecy_diff}
    f.close()


g = open(f"{work_dir}/variant.TF_motif.overlap.combined.motif_disruption.match_threshold.uniq","w")
g.write(header)
for pair in pair2info:
    g.write("\t".join(pair2info[pair]["info"]) + "\n")
g.close()