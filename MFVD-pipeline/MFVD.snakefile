#----set the path of configure file----#
configfile_dir = "."
configfile: configfile_dir + "/config.yaml"

#----set rule all
rule all:
    input:
        config["projectDir"].rstrip("/") + "/tf_regulation/tfieQTL_result_merged/tfieQTL.all.top_eVariant_for_each_TF_target_pair.withGeneName"
        config["projectDir"].rstrip("/") + "/ASoC_result/all.txt"

#----ATAC-seq data alignment and peak calling----#
rule bowtie2_alignment:
    input:
        lambda wildcards: config["fastqFiles"][wildcards.sample + "_1"],
        lambda wildcards: config["fastqFiles"][wildcards.sample + "_2"]
    output:
        temp(config["projectDir"].rstrip("/") + "/alignment_result/sam/{sample}.sam")
    threads: config["Threads"]["bowtie2"]
    shell:
        f'{config["softwarePath"]["bowtie2"]} -p {{threads}} -x {config["bowtie2IndexPrefix"]} -1 {{input[0]}} -2 {{input[1]}} -S {{output}}'
 
rule sam2bam_and_filter_bam:
    input:
        config["projectDir"].rstrip("/") + '/alignment_result/sam/{sample}.sam'
    output:
        temp(config["projectDir"].rstrip("/") + '/alignment_result/filtered_bam/{sample}_filtered.bam')
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} view -@ {{threads}} -b -q 30 {{input}} -o {{output}}'

rule sort_bam_files:
    input:
        config["projectDir"].rstrip("/") + "/alignment_result/filtered_bam/{sample}_filtered.bam"
    output:
        temp(config["projectDir"].rstrip("/") + "/alignment_result/sorted_bam/{sample}_filtered_sorted.bam")
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} sort -@ {{threads}} -o {{output}} {{input}}'

rule remove_duplication:
    input:
        config["projectDir"].rstrip("/") + "/alignment_result/sorted_bam/{sample}_filtered_sorted.bam"
    output:
        config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam/{sample}_filtered_sorted_rmdup.bam",
        temp(config["projectDir"].rstrip("/") + '/alignment_result/{sample}_marked_dup_metrics.txt')
    shell:
        f'{config["softwarePath"]["picard"]} MarkDuplicates I={{input}} O={{output[0]}} M={{output[1]}} REMOVE_DUPLICATES=True'

rule bam_index:
    input:
        config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam/{sample}_filtered_sorted_rmdup.bam"
    output:
        config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam/{sample}_filtered_sorted_rmdup.bam.bai"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} index -@ {{threads}} {{input}}'

rule macs2_call_peak_shift_model:
    input:
        config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam/{sample}_filtered_sorted_rmdup.bam"
    output:
        # config["projectDir"].rstrip("/") + "/peak_results/{sample,r'^(?!.*all\_samples).*$'}_peaks.narrowPeak"
        config["projectDir"].rstrip("/") + "/peak_results/{sample}_peaks.narrowPeak"
    wildcard_constraints:
        # sample = r'^(?!.*all).*$'
        sample = '.*ATAC.*'
        # config["regularExpressionOfSampleNames"]
    params:
        config["projectDir"].rstrip("/") + "/peak_results/{sample}"
    shell:
        f'{config["softwarePath"]["MACS2"]} callpeak --nomodel -t {{input}} -g hs -n {{params}} -q 0.05 --extsize 200 --shift -100 --keep-dup all'

rule extract_ATAC_peak_present_in_multi_sample:
    input:
        expand(config["projectDir"].rstrip("/") + "/peak_results/{sample}_peaks.narrowPeak",sample=config["samples"])
    output:
        config["projectDir"].rstrip("/") + "/peak_results_summary/ATAC_peak_in_multiple_samples.merged.bed"
    params:
        sampleCount = config["sampleCount"],
        inputDir = config["projectDir"].rstrip("/") + "/peak_results",
        outputDir = config["projectDir"].rstrip("/") + "/peak_results_summary",
        bedopsPath = config["softwarePath"]["bedops"],
        bedtoolsPath = config["softwarePath"]["bedtools"]
    shell:
        'bash ./scripts/extract_ATAC_peak_present_in_multi_sample.sh {params.bedopsPath} {params.bedtoolsPath} {params.sampleCount} {params.inputDir} {params.outputDir} {input}'

rule merge_all_sample_bams:
    input:
        bam = expand(config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam/{sample}_filtered_sorted_rmdup.bam",sample=config["samples"]),
        index = expand(config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam/{sample}_filtered_sorted_rmdup.bam.bai",sample=config["samples"])
    output:
        config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam_merged/all_filtered_sorted_rmdup.bam"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} merge -@ {{threads}} -f {{output}} {{input.bam}}'

rule merged_bam_index:
    input:
        config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam_merged/all_filtered_sorted_rmdup.bam"
    output:
        config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam_merged/all_filtered_sorted_rmdup.bam.bai"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} index -@ {{threads}} {{input}}'

rule macs2_call_peak_shift_model_for_merged_bam:
    input:
        config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam_merged/all_filtered_sorted_rmdup.bam"
    output:
        config["projectDir"].rstrip("/") + "/peak_results/all_samples_peaks.narrowPeak"
    params:
        config["projectDir"].rstrip("/") + "/peak_results/all_samples"
    shell:
        f'{config["softwarePath"]["MACS2"]} callpeak --nomodel -t {{input}} -g hs -n {{params}} -q 0.05 --extsize 200 --shift -100 --keep-dup all --call-summits'

rule extract_final_peaks:
    input:
        config["projectDir"].rstrip("/") + "/peak_results/all_samples_peaks.narrowPeak",
        config["projectDir"].rstrip("/") + "/peak_results_summary/ATAC_peak_in_multiple_samples.merged.bed"
    output:
        config["projectDir"].rstrip("/") + "/peak_results_summary/all_samples_peaks_in_mutiple_samples.narrowPeak"
    shell:
        f'{config["softwarePath"]["bedtools"]} intersect -u -a {{input[0]}} -b {{input[1]}} > {{output}}'

rule extract_uniq_peak_region:
    input:
        config["projectDir"].rstrip("/") + "/peak_results_summary/all_samples_peaks_in_mutiple_samples.narrowPeak"
    output:
        config["projectDir"].rstrip("/") + "/peak_results_summary/all_samples_peaks_in_mutiple_samples.narrowPeak.uniq"
    shell:
        'cut -f 1-3 {input} | sort | uniq > {output}'

#----rules to produce input files for TF footprint analysis----#
rule extract_motifs:
    input:
        jasparMotif = config["motifFile"]["jaspar"],
        hocomocoMotif = config["motifFile"]["HOCOMOCO"],
        expressedGene = config["geneExpressionFile"],
    output:
        jaspar = config["projectDir"].rstrip("/") + "/motifFile/jaspar_motifs_of_expressed_TF.txt",
        homocomo = config["projectDir"].rstrip("/") + "/motifFile/hocomoco_motifs_of_expressed_TF.txt"
    params: config["expressionThreshold"]
    script:
        "./scripts/extract_motifs.py"

rule produce_snp_file_as_SNPsplit_input:
    input:
        bimFile = config["bimFile"],
        ref = config["referenceGenomeFastaFile"]
    output:
        expand("SNPs_from_user/chr{chr_}.txt",chr_ = list(range(1,23)))
    params:
        "SNPs_from_user"
    threads:
        config["Threads"]["produceSnpFiles"]
    script:
        "./scripts/produce_snp_files_as_SNPsplit_input.py"

rule produce_genome_sequence_with_snps_trasfermed:
    input:
        expand("SNPs_from_user/chr{chr_}.txt",chr_ = list(range(1,23)))
    output:
        expand("from_user_full_sequence/chr{chr_}.SNPs_introduced.fa",chr_ = list(range(1,23)) + ["X","Y","M"])
    shell:
        f'{config["softwarePath"]["SNPsplit_genome_preparation"]} --strain from_user '
        f'--reference_genome {config["referenceGenomeFastaFileDir"]} '
        f'--skip_filtering --full_sequence {config["referenceGenomeFastaFilePrefix"]}'

rule merge_variants_transfered_reference_genome:
    input:
        expand("from_user_full_sequence/chr{chr_}.SNPs_introduced.fa",chr_ = list(range(1,23)) + ["X","Y","M"])
    output:
        f'from_user_full_sequence/{config["referenceGenomeFastaFilePrefix"]}.SNPs_introduced.fa'
    params:
        "from_user_full_sequence"
    script:
        "./scripts/merge_variants_transfered_reference_genome.py"

#----rules for TF footprint analysis----#
rule run_TOBAIA_ATACorrect:
    input:
        ref = config["referenceGenomeFastaFile"],
        bam = config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam_merged/all_filtered_sorted_rmdup.bam",
        bamindex = config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam_merged/all_filtered_sorted_rmdup.bam.bai",
        peaks = config["projectDir"].rstrip("/") + "/peak_results_summary/all_samples_peaks_in_mutiple_samples.narrowPeak.uniq",
        blackList = config["genomeBlackList"]
    output:
        config["projectDir"].rstrip("/") + "/tobias_footprint/01_ATACorrect/all_samples_corrected.bw"
    params:
        config["projectDir"].rstrip("/") + "/tobias_footprint/01_ATACorrect"
    threads: config["Threads"]["TOBIAS"]
    shell:
        f'{config["softwarePath"]["TOBIAS"]} ATACorrect --bam {{input.bam}} --genome {{input.ref}} '
        f'--peaks {{input.peaks}} --blacklist {{input.blackList}} --outdir {{params}} --prefix all_samples --cores {{threads}}'

rule run_TOBAIA_FootprintScores:
    input:
        bw = config["projectDir"].rstrip("/") + "/tobias_footprint/01_ATACorrect/all_samples_corrected.bw",
        peaks = config["projectDir"].rstrip("/") + "/peak_results_summary/all_samples_peaks_in_mutiple_samples.narrowPeak.uniq"
    output:
        config["projectDir"].rstrip("/") + "/tobias_footprint/02_FootprintScores/all_samples_footprints.bw"
    threads: config["Threads"]["TOBIAS"]
    shell:
        f'{config["softwarePath"]["TOBIAS"]} FootprintScores --signal {{input.bw}} --regions {{input.peaks}} --output {{output}} --cores {{threads}}'

rule run_TOBAIA_BINDetect:
    input:
        bw = config["projectDir"].rstrip("/") + "/tobias_footprint/02_FootprintScores/all_samples_footprints.bw",
        ref = config["referenceGenomeFastaFile"],
        alter_ref = f'from_user_full_sequence/{config["referenceGenomeFastaFilePrefix"]}.SNPs_introduced.fa',
        peaks = config["projectDir"].rstrip("/") + "/peak_results_summary/all_samples_peaks_in_mutiple_samples.narrowPeak.uniq",
        jaspar = config["projectDir"].rstrip("/") + "/motifFile/jaspar_motifs_of_expressed_TF.txt",
        homocomo = config["projectDir"].rstrip("/") + "/motifFile/hocomoco_motifs_of_expressed_TF.txt"
    output:
        directory(config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_jaspar_originalRef"),
        directory(config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_jaspar_snpsIntroducedRef"),
        directory(config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_HOCOMOCO_originalRef"),
        directory(config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_HOCOMOCO_snpsIntroducedRef")
    params:
        outputDir = config["projectDir"].rstrip("/") + "/tobias_footprint",
        tobiasPath = config["softwarePath"]["TOBIAS"]
    threads: config["Threads"]["TOBIAS"]
    shell:
        "bash ./scripts/run_TOBAIA_BINDetect.sh {params.tobiasPath} {threads} {input.bw} {input.ref} {input.alter_ref} {input.peaks} {input.jaspar} {input.homocomo} {params.outputDir}"

rule merge_footprints_for_each_database:
    input:
        config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_jaspar_originalRef",
        config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_jaspar_snpsIntroducedRef",
        config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_HOCOMOCO_originalRef",
        config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_HOCOMOCO_snpsIntroducedRef"
    output:
        config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_combined.new_footprint_count.summary",
        directory(config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_jaspar_combined"),
        directory(config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_HOCOMOCO_combined")
    params:
        config["projectDir"].rstrip("/") + "/tobias_footprint"
    script:
        "./scripts/merge_footprint_results_from_original_and_snpsintroducedRef.py"

rule merge_footprint_to_one_bed_file:
    input:
        config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_jaspar_combined",
        config["projectDir"].rstrip("/") + "/tobias_footprint/03_BINDetect_HOCOMOCO_combined",
        config["geneName2id_anno_file"]
    output:
        config["projectDir"].rstrip("/") + "/tobias_footprint/04_all_motifs_in_one_file/jaspar.all_motif_footprint.bed",
        config["projectDir"].rstrip("/") + "/tobias_footprint/04_all_motifs_in_one_file/HOCOMOCO.all_motif_footprint.bed"
    params:
        inputDir = config["projectDir"].rstrip("/") + "/tobias_footprint",
        outputDir = config["projectDir"].rstrip("/") + "/tobias_footprint/04_all_motifs_in_one_file"
    script:
        "./scripts/merge_footprint_of_all_TFs_to_one_file.py"

#----rules for TF binding disruption analysis----#
rule get_variants_footprint_overlap:
    input:
        config["variantsBedFile"],
        config["projectDir"].rstrip("/") + "/tobias_footprint/04_all_motifs_in_one_file/jaspar.all_motif_footprint.bed",
        config["projectDir"].rstrip("/") + "/tobias_footprint/04_all_motifs_in_one_file/HOCOMOCO.all_motif_footprint.bed"
    output:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.jaspar",
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.HOCOMOCO"
    params:
        inputDir = config["projectDir"].rstrip("/") + "/tobias_footprint/04_all_motifs_in_one_file",
        outputDir = config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap",
        bedtools = config["softwarePath"]["bedtools"]
    shell:
        "bash ./scripts/get_variants_footprint_overlap.sh {params.bedtools} {input[0]} {params.inputDir} {params.outputDir}"

rule get_motif_infromation_content_of_and_frequeny_diff:
    input:
        jaspar_IC = config["motifFile"]["jaspar_IC"],
        homocomo_IC = config["motifFile"]["HOCOMOCO_IC"],
        jaspar_overlap = config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.jaspar",
        hocomoco_overlap = config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.HOCOMOCO"     
    output:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.jaspar.motif_disruption",
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.HOCOMOCO.motif_disruption"    
    params:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap"
    script:
        "./scripts/get_motif_infromation_content_of_and_frequeny_diff.py"

rule extract_eVariants_with_motif_disruption_above_threshold:
    input:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.jaspar.motif_disruption",
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.HOCOMOCO.motif_disruption"
    output:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.jaspar.motif_disruption.match_threshold",
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.HOCOMOCO.motif_disruption.match_threshold"
    params:
        fileDir = config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap",
        ic_threshold = config["motifDisruptionThreshold"]["information_content_of_position"],
        frequency_threshold = config["motifDisruptionThreshold"]["frequency_diff_of_two_allele"]
    shell:
        "bash ./scripts/extract_eVariants_with_motif_disruption_above_threshold.sh {params.fileDir} {params.ic_threshold} {params.frequency_threshold}"

rule retain_the_top_frequency_diff_for_each_TF_target_variant_pair:
    input:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.jaspar.motif_disruption.match_threshold",
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_footprint.overlap.HOCOMOCO.motif_disruption.match_threshold"
    output:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_motif.overlap.combined.motif_disruption.match_threshold.uniq"
    params:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap"
    script:
        "./scripts/retain_the_top_frequency_diff_for_each_TF_target_variant_pair.py"

rule prepare_input_for_tfieQTL:
    input:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/variant.TF_motif.overlap.combined.motif_disruption.match_threshold.uniq"
    output:
        expand(config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/splited_evariant_motif_overlap_files/variant.TF_motif.overlap_{batch}",batch=list(range(1,11)))
    params:
        config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/splited_evariant_motif_overlap_files"
    script:
        "./scripts/prepare_input_for_tfieQTL.py"

rule tfieQTL:
    input:
        eVariant2TF = config["projectDir"].rstrip("/") + "/tf_regulation/variant_footprint_overlap/splited_evariant_motif_overlap_files/variant.TF_motif.overlap_{batch}",
        exp = config["expressionBedFile"],
        cov = config["covariateFile"],
        samplelist = config["sampleList"]
    output:
        config["projectDir"].rstrip("/") + "/tf_regulation/tfieQTL_result/tfieQTL.{batch}"
    conda:
        "cieQTL"
    script:
        "./scripts/tfieQTL_for_snakemake.R"

rule merge_tfieQTL_result:
    input:
        expand(config["projectDir"].rstrip("/") + "/tf_regulation/tfieQTL_result/tfieQTL.{batch}",batch=list(range(1,11)))
    output:
        config["projectDir"].rstrip("/") + "/tf_regulation/tfieQTL_result_merged/tfieQTL.all"
    shell:
        "cat {input} | sed '/interactP/d' | sort | uniq | sort -gk10,10 "
        "| sed '1i TF\\ttaget_gene\\tSNP\\tAICDiff\\tchiSquarePvalue\\tEstimate\\tSE\\tdf\\tt_value\\tinteractP' "
        "> {output}"

rule filter_merged_tfieQTL_and_retain_the_top_variants_for_each_TF_target_pair:
    input:
        config["projectDir"].rstrip("/") + "/tf_regulation/tfieQTL_result_merged/tfieQTL.all"
    output:
        config["projectDir"].rstrip("/") + "/tf_regulation/tfieQTL_result_merged/tfieQTL.all.top_eVariant_for_each_TF_target_pair"
    script:
        "./scripts/filter_merged_tfieQTL_and_retain_the_top_variants_for_each_TF_target_pair.py"

rule annotate_geneID:
    input:
        config["geneName2id_anno_file"],
        config["projectDir"].rstrip("/") + "/tf_regulation/tfieQTL_result_merged/tfieQTL.all.top_eVariant_for_each_TF_target_pair"
    output:
        config["projectDir"].rstrip("/") + "/tf_regulation/tfieQTL_result_merged/tfieQTL.all.top_eVariant_for_each_TF_target_pair.withGeneName"
    script:
        "./scripts/annotate_geneID.py"

#----rules for allele-specific open chromotin SNPs analysis ----#
rule transfer_vcf_files_to_h5_files:
    input:
        vcf = expand(config["vcfFilesDir"] + "/chr{chr_}.vcf.gz",chr_ = list(range(1,23))),
        chr_length = config["chromosomeLengthFile"]
    output:
        haplotype = config["projectDir"].rstrip("/") + "/WASP_correct/h5_files/haplotypes.h5",
        snp_index = config["projectDir"].rstrip("/")+ "/WASP_correct/h5_files/snp_index.h5",
        snp_tab = config["projectDir"].rstrip("/") + "/WASP_correct/h5_files/snp_tab.h5"
    shell:
        f'{config["softwarePath"]["snp2h5"]} --chrom {{input.chr_length}} --format vcf --haplotype {{output.haplotype}} --snp_index {{output.snp_index}} --snp_tab {{output.snp_tab}} {{input.vcf}}'

rule find_intersecting_snps:
    input:
        bam = config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam/{sample}_filtered_sorted_rmdup.bam",
        index = config["projectDir"].rstrip("/") + "/alignment_result/rmdup_bam/{sample}_filtered_sorted_rmdup.bam.bai",
        haplotype = config["projectDir"].rstrip("/") + "/WASP_correct/h5_files/haplotypes.h5",
        snp_index = config["projectDir"].rstrip("/") + "/WASP_correct/h5_files/snp_index.h5",
        snp_tab = config["projectDir"] + "/WASP_correct/h5_files/snp_tab.h5"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.keep.bam",
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.to.remap.bam",
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.remap.fq1.gz",
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.remap.fq2.gz"
    params:
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps"
    shell:
        f'python {config["softwarePath"]["find_intersecting_snps_modified"]} --is_paired_end --is_sorted '
        f'--output_dir {{params}} --haplotype {{input.haplotype}} --snp_index {{input.snp_index}} --snp_tab {{input.snp_tab}} {{input.bam}}'

rule transfer_fake_bam_to_real_bam:
    input:
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.keep.bam"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.keep.real.bam"
    threads:
        config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} view -@ {{threads}} -b -o {{output}} {{input}}'

rule bowtie2_remap:
    input:
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.remap.fq1.gz",
        config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.remap.fq2.gz"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/rmake_bam/{sample}.remap.sam"
    shell:
        f'{config["softwarePath"]["bowtie2"]} -p {{threads}} -x {config["bowtie2IndexPrefix"]} -1 {{input[0]}} -2 {{input[1]}} -S {{output}}'

rule sort_remap_sam_file:
    input:
        config["projectDir"].rstrip("/") + "/WASP_correct/rmake_bam/{sample}.remap.sam"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/rmake_bam/{sample}.remap.sorted.bam"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} sort -@ {{threads}} -o {{output}} {{input}}'

rule remap_bam_index:
    input:
        config["projectDir"].rstrip("/") + "/WASP_correct/rmake_bam/{sample}.remap.sorted.bam"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/rmake_bam/{sample}.remap.sorted.bam.bai"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} index -@ {{threads}} {{input}}'

rule filter_remapped_reads:
    input:
        to_remap_bam = config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.to.remap.bam",
        remap_bam = config["projectDir"].rstrip("/") + "/WASP_correct/rmake_bam/{sample}.remap.sorted.bam"
    output:
        keep_bam = config["projectDir"].rstrip("/") + "/WASP_correct/filter_remapped_reads/{sample}.remap.sorted.keep.bam"
    shell:
        f'python {config["softwarePath"]["filter_remapped_reads"]} {{input.to_remap_bam}} {{input.remap_bam}} {{output.keep_bam}}'

rule merge_wasp_corrected_bams:
    input:
        first_bam = config["projectDir"].rstrip("/") + "/WASP_correct/find_intersecting_snps/{sample}_filtered_sorted_rmdup.keep.real.bam",
        second_bam = config["projectDir"].rstrip("/") + "/WASP_correct/filter_remapped_reads/{sample}.remap.sorted.keep.bam"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.bam"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} merge -@ {{threads}} -f {{output}} {{input.first_bam}} {{input.second_bam}}'

rule sort_merged_wasp_corrected_bam:
    input:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.bam"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.bam"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} sort -@ {{threads}} -o {{output}} {{input}}'

rule mrged_wasp_bam_index:
    input:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.bam"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.bam.bai"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} index -@ {{threads}} {{input}}'

rule add_RG_filed_to_bam_file:
    input:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.bam",
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.bam.bai"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.withRG.bam"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} addreplacerg -@ {{threads}} -r "@RG\\tID:{{wildcards.sample}}\\tPL:ILLUMINA\\tLB:{{wildcards.sample}}\\tSM:{{wildcards.sample}}" -o {{output}} {{input}}'

rule index_RG_contained_bam:
    input:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.withRG.bam"
    output:
        config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.withRG.bam.bai"
    threads: config["Threads"]["samtools"]
    shell:
        f'{config["softwarePath"]["samtools"]} index -@ {{threads}} {{input}}'

rule run_gatk_HaplotypeCaller:
    input:
        bam = config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.withRG.bam",
        index = config["projectDir"].rstrip("/") + "/WASP_correct/merged_bam/{sample}.all.keep.sorted.withRG.bam.bai",
        ref = config["gatk_resource"]["gatk_ref_genome"],
        dbsnp = config["gatk_resource"]["dbsnp_vcf_file"]
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/GVCF_results_from_HaplotypeCaller/{sample}.gvcf.gz"
    shell:
        f'{config["softwarePath"]["gatk"]} --java-options "-Xmx4g" HaplotypeCaller -R {{input.ref}} --dbsnp {{input.dbsnp}} -I {{input.bam}} -O {{output}} -ERC GVCF'

rule genotype_gvcf:
    input:
        gvcf = config["projectDir"].rstrip("/") + "/GATK_result/GVCF_results_from_HaplotypeCaller/{sample}.gvcf.gz",
        ref = config["gatk_resource"]["gatk_ref_genome"]
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/raw_VCF_result/{sample}.raw_variant.vcf.gz"
    shell:
        f'{config["softwarePath"]["gatk"]} --java-options "-Xmx4g" GenotypeGVCFs -R {{input.ref}} -V {{input.gvcf}} -O {{output}}'

rule Build_recalibration_model_for_InDel:
    input:
        vcf = config["projectDir"].rstrip("/") + "/GATK_result/raw_VCF_result/{sample}.raw_variant.vcf.gz",
        ref = config["gatk_resource"]["gatk_ref_genome"],
        gold_standard_Indel = config["gatk_resource"]["gold_standard_Indel"]
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_InDel.recal",
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_InDel.tranches",
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_InDel.plots.R"
    conda:
        "cieQTL"
    shell:
        f'{config["softwarePath"]["gatk"]} VariantRecalibrator --trust-all-polymorphic '
        f'-max-gaussians 4 -R {{input.ref}} -V {{input.vcf}} '
        f'--resource:mills,known=true,training=true,truth=true,prior=12.0 {{input.gold_standard_Indel}} '
        '-an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR '
        '-mode INDEL '
        '-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 '
        f'-O {{output[0]}} --tranches-file {{output[1]}} --rscript-file {{output[2]}}'

rule Build_recalibration_model_for_SNP:
    input:
        vcf = config["projectDir"].rstrip("/") + "/GATK_result/raw_VCF_result/{sample}.raw_variant.vcf.gz",
        ref = config["gatk_resource"]["gatk_ref_genome"],
        hampam_known = config["gatk_resource"]["hampam_known"],
        genome_1000_omni = config["gatk_resource"]["genome_1000_omni"],
        high_confidence = config["gatk_resource"]["high_confidence"],
        db_snp_known = config["gatk_resource"]["dbsnp_vcf_file"]
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_SNP.recal",
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_SNP.tranches",
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_SNP.plots.R"
    conda:
        "cieQTL"
    shell:
        f'{config["softwarePath"]["gatk"]} VariantRecalibrator --trust-all-polymorphic '
        f'--max-gaussians 4 -R {{input.ref}} -V {{input.vcf}} '
        f'--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {{input.hampam_known}} '
        f'--resource:omni,known=false,training=true,truth=true,prior=12.0 {{input.genome_1000_omni}} '
        f'--resource:1000G,known=false,training=true,truth=false,prior=10.0 {{input.high_confidence}} '
        f'--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {{input.db_snp_known}} '
        '-an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR '
        '-mode SNP '
        '-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 '
        f'-O {{output[0]}} --tranches-file {{output[1]}} --rscript-file {{output[2]}}'

rule ApplyVQSR_for_InDel:
    input:
        config["projectDir"].rstrip("/") + "/GATK_result/raw_VCF_result/{sample}.raw_variant.vcf.gz",
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_InDel.recal",
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_InDel.tranches",
        config["gatk_resource"]["gatk_ref_genome"]
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_recalibrated.vcf.gz"
    shell:
        f'{config["softwarePath"]["gatk"]} ApplyVQSR -mode INDEL '
        f'-R {{input[3]}} '
        f'-V {{input[0]}} --recal-file {{input[1]}} --tranches-file {{input[2]}} '
        '--truth-sensitivity-filter-level 99.0 '
        '--create-output-variant-index true '
        f'-O {{output}}'

rule ApplyVQSR_for_SNP:
    input:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_recalibrated.vcf.gz",
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_SNP.recal",
        config["projectDir"].rstrip("/") + "/GATK_result/VQSR_results/{sample}_SNP.tranches",
        config["gatk_resource"]["gatk_ref_genome"]
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_SNP_recalibrated.vcf.gz"
    shell:
        f'{config["softwarePath"]["gatk"]} ApplyVQSR -mode SNP '
        f'-R {{input[3]}} '
        f'-V {{input[0]}} --recal-file {{input[1]}} --tranches-file {{input[2]}} '
        '--truth-sensitivity-filter-level 99.0 '
        '--create-output-variant-index true '
        f'-O {{output}}'

rule remove_Indel_from_vcf:
    input:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_SNP_recalibrated.vcf.gz"
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_SNP_recalibrated.PASS_only.vcf.gz"
    shell:
        f'{config["softwarePath"]["vcftools"]} --gzvcf {{input}} --remove-filtered-all --remove-indels --max-alleles 2 --recode --stdout | gzip -c > {{output}}'

rule remove_snp_not_in_peak_regions:
    input:
        vcf = config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_SNP_recalibrated.PASS_only.vcf.gz",
        peak = config["projectDir"].rstrip("/") + "/peak_results_summary/all_samples_peaks_in_mutiple_samples.narrowPeak.uniq"
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_SNP_recalibrated.PASS_only.within_peak.vcf.gz"
    shell:
        f'{config["softwarePath"]["vcftools"]} --gzvcf {{input.vcf}} --bed {{input.peak}} --recode --stdout | gzip -c > {{output}}'

rule extract_heterozygote:
    input:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_SNP_recalibrated.PASS_only.within_peak.vcf.gz"
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_SNP_recalibrated.PASS_only.within_peak.heterozygosis.vcf.gz"
    params:
        config["softwarePath"]["tabix"]
    script:
        "./scripts/extract_heterozygote.py"

rule merge_filtered_vcf_files:
    input:
        expand(config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/{sample}_InDel_SNP_recalibrated.PASS_only.within_peak.heterozygosis.vcf.gz",sample=config["samples"])
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/all_samples_InDel_SNP_recalibrated.PASS_only.within_peak.heterozygosis.vcf.gz"
    shell:
        f'{config["softwarePath"]["bcftools"]} merge {{input}} | gzip -c > {{output}}'

rule allele_effect_direction_check:
    input:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/all_samples_InDel_SNP_recalibrated.PASS_only.within_peak.heterozygosis.vcf.gz"
    output:
        config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/all_samples_InDel_SNP_recalibrated.PASS_only.within_peak.heterozygosis.vcf.summary"
    script:
        "./scripts/allele_effect_direction_check.py"

rule binom_test:
    input:
        vcf_summary = config["projectDir"].rstrip("/") + "/GATK_result/filetered_VCF_results/all_samples_InDel_SNP_recalibrated.PASS_only.within_peak.heterozygosis.vcf.summary"
    output:
        all_ = config["projectDir"].rstrip("/") + "/ASoC_result/all.txt",
        p = config["projectDir"].rstrip("/") + "/ASoC_result/significant_pvalue_0.05.txt",
        fdr_1 = config["projectDir"].rstrip("/") + "/ASoC_result/significant_fdr_0.1.txt",
        fdr_2 = config["projectDir"].rstrip("/") + "/ASoC_result/significant_fdr_0.05.txt"
    conda:
        "vcftools"
    script:
        "./scripts/binom_test.R"