import os

tobias_dir = snakemake.params[0]

new_footprint_count_summary = open(f"{tobias_dir}/03_BINDetect_combined.new_footprint_count.summary","w")
new_footprint_count_summary.write("motifDatabase\tmotif\tnew_footprint_count_from_snps_introduced_ref\tratio_of_new_footprint_count_compared_to_that_from_original_ref\n")

for motifDatabase in ["HOCOMOCO","jaspar"]:
    outdir = f"{tobias_dir}/03_BINDetect_{motifDatabase}_combined"
    os.system(f"mkdir -p {outdir}")
    for motif in os.listdir(f"{tobias_dir}/03_BINDetect_{motifDatabase}_originalRef"):
        if not motif.endswith("txt") and not motif.endswith("pdf") and not motif.endswith("xlsx") and "snakemake" not in motif:
            # write footprint identified from original genome
            footprint_from_originalRef = set()
            g = open(f"{outdir}/{motif}_all_samples_footprints_bound.bed","w")
            f = open(f"{tobias_dir}/03_BINDetect_{motifDatabase}_originalRef/{motif}/beds/{motif}_all_samples_footprints_bound.bed")
            for line in f:
                g.write(line.strip() + "\tfrom_originalRef\n")
                lines = line.strip().split("\t")
                footprint_from_originalRef.add("_".join(lines[:3]))
            f.close()

            # write footprint identified only from snps introduced genome
            new_footprint_count = 0
            f = open(f"{tobias_dir}/03_BINDetect_{motifDatabase}_snpsIntroducedRef/{motif}/beds/{motif}_all_samples_footprints_bound.bed")
            for line in f:
                lines = line.strip().split("\t")
                if "_".join(lines[:3]) not in footprint_from_originalRef:
                    g.write(line.strip() + "\tfrom_snpsIntroducedRef\n")
                    new_footprint_count += 1
            f.close()
            g.close()

            new_footprint_count_summary.write("\t".join([motifDatabase,motif,str(new_footprint_count),str(new_footprint_count/len(footprint_from_originalRef))]) + "\n")

new_footprint_count_summary.close()
