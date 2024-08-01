import os

#----define function to load motif information content----#
def loadMotif(tfDatabase):
    motif2info = {}

    if tfDatabase == "ENCODE":
        f = open("/home/TW/data/TFs/ENCODE_2013/ENCODE_2013.txt.information_content")
        for line in f:
            if line.startswith(">"):
                motif = line.rstrip("\n").lstrip(">").split()[0]
                motif2info[motif] = {}
            else:
                index,ic,A,C,G,T = line.strip().split()
                motif2info[motif][index] = {"ic":float(ic),"A":float(A),"C":float(C),"G":float(G),"T":float(T)}
        f.close()
    
    elif tfDatabase == "jaspar":
        f = open("/home/TW/data/TFs/jaspar2022/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt.information_content")
        for line in f:
            if line.startswith(">"):
                motif = line.rstrip("\n").split()[0][1:]
                motif2info[motif] = {}
            else:
                index,ic,A,C,G,T = line.strip().split()
                motif2info[motif][index] = {"ic":float(ic),"A":float(A),"C":float(C),"G":float(G),"T":float(T)}
        f.close()

    elif tfDatabase == "HOCOMOCO":
        f = open("/home/TW/data/TFs/HOCOMOCO_v11/HOCOMOCOv11_core_HUMAN_mono_jaspar_format.tfName.txt.information_content")
        for line in f:
            if line.startswith(">"):
                motif = line.rstrip("\n").split()[0][1:]
                motif2info[motif] = {}
            else:
                index,ic,A,C,G,T = line.strip().split()
                motif2info[motif][index] = {"ic":float(ic),"A":float(A),"C":float(C),"G":float(G),"T":float(T)}
    
    return motif2info


work_dir = snakemake.params[0]
letterTransfer = {"A":"T", "T":"A", "C":"G", "G":"C"}

#----extract infromation content of tfieQTL significant variant----#
for tfDatabase in ["jaspar","HOCOMOCO"]:
    motif2info = loadMotif(tfDatabase)
    f = open(f"{work_dir}/eVariant.TF_footprint.overlap.{tfDatabase}")
    g = open(f"{work_dir}/eVariant.TF_footprint.overlap.{tfDatabase}.motif_disruption","w")
    g.write(f.readline().rstrip("\n") + "\t" +  "\t".join(["motif_informat_content","major_frequency","minor_frequency","frequecy_diff"]) + "\n")
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
            
            if motif_direction == "+":
                disrupted_motid_index = variant_pos-motif_start+1
                letter_infomation_content = motif2info[motif][str(disrupted_motid_index)]["ic"]
                frequency_allele1 = motif2info[motif][str(disrupted_motid_index)][allele1]
                frequency_allele2 = motif2info[motif][str(disrupted_motid_index)][allele2]
                major_frequency = max([frequency_allele1,frequency_allele2])
                minor_frequency = min([frequency_allele1,frequency_allele2])
                frequency_diff = abs(frequency_allele1 - frequency_allele2)


            if motif_direction == "-":
                disrupted_motid_index = motif_end-variant_pos+1
                letter_infomation_content = motif2info[motif][str(disrupted_motid_index)]["ic"]
                frequency_allele1 = motif2info[motif][str(disrupted_motid_index)][letterTransfer[allele1]]
                frequency_allele2 = motif2info[motif][str(disrupted_motid_index)][letterTransfer[allele2]]
                major_frequency = max([frequency_allele1,frequency_allele2])
                minor_frequency = min([frequency_allele1,frequency_allele2])
                frequency_diff = abs(frequency_allele1 - frequency_allele2)
            
            g.write(line.rstrip("\n") + "\t" + "\t".join([str(letter_infomation_content),str(major_frequency),str(minor_frequency),str(frequency_diff)]) + "\n")
    f.close()
    g.close()
