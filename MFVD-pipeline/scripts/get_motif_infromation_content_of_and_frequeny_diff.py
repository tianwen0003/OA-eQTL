import os

#----define function to load motif information content----#
def loadMotif(tfDatabase):
    motif2info = {}

    if tfDatabase == "jaspar":
        f = open(snakemake.input[0])
        for line in f:
            if line.startswith(">"):
                motif = line.rstrip("\n").split()[0][1:]
                motif2info[motif] = {}
            else:
                index,ic,A,C,G,T = line.strip().split()
                motif2info[motif][index] = {"ic":float(ic),"A":float(A),"C":float(C),"G":float(G),"T":float(T)}
        f.close()

    elif tfDatabase == "HOCOMOCO":
        f = open(snakemake.input[1])
        for line in f:
            if line.startswith(">"):
                motif = line.rstrip("\n").split()[0][1:]
                motif2info[motif] = {}
            else:
                index,ic,A,C,G,T = line.strip().split()
                motif2info[motif][index] = {"ic":float(ic),"A":float(A),"C":float(C),"G":float(G),"T":float(T)}
    
    return motif2info


work_dir = snakemake.params[0]

#----extract infromation content of tfieQTL significant variant----#
for tfDatabase in ["jaspar","HOCOMOCO"]:
    motif2info = loadMotif(tfDatabase)
    f = open(f"{work_dir}/variant.TF_footprint.overlap.{tfDatabase}")
    g = open(f"{work_dir}/variant.TF_footprint.overlap.{tfDatabase}.motif_disruption","w")
    g.write(f.readline().rstrip("\n") + "\t" +  "\t".join(["motif_informat_content","major_allele_in_motif","minor_allele_in_motif","major_frequency","minor_frequency","frequecy_diff"]) + "\n")
    for line in f:
        lines = line.strip().split("\t")
        
        variant_pos = int(lines[1])
        allele1,allele2 = lines[6:8]

        motif = "_".join(lines[13].split("_")[1:])
        if len(allele1) ==1 and len(allele2) == 1:
            motif_chr = lines[10]
            motif_start = int(lines[11])
            motif_end = int(lines[12])

            motif_position = "_".join(lines[10:13])

            motif_direction = lines[16]

            disrupted_motid_index = variant_pos-motif_start+1

            if motif_direction == "-":
                disrupted_motid_index = motif_end-variant_pos+1

            letter_infomation_content = motif2info[motif][str(disrupted_motid_index)]["ic"]
            frequency_allele1 = motif2info[motif][str(disrupted_motid_index)][allele1]
            frequency_allele2 = motif2info[motif][str(disrupted_motid_index)][allele2]
            major_frequency = max([frequency_allele1,frequency_allele2])
            minor_frequency = min([frequency_allele1,frequency_allele2])
            frequency_diff = abs(frequency_allele1 - frequency_allele2)

            major_allele = allele1
            minor_allele = allele2

            if frequency_allele1 < frequency_allele2:
                major_allele = allele2
                minor_allele = allele1
            
            g.write(line.rstrip("\n") + "\t" + "\t".join([str(letter_infomation_content),major_allele,minor_allele,str(major_frequency),str(minor_frequency),str(frequency_diff)]) + "\n")
    f.close()
    g.close()