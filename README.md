# Introduce
This repository contains the codes of our established Multi-level functional variants-deciphering pipeline based on snakemake, which was used for functional fine-mapping of eQTL SNPs our OA-eQTL project. 

The pipeline is organized as a snakemake workflow in which invoked some custom python and R scripts, so we assume the users are familiar with the basic principles of snakemake, python and R.

# How to use
## Create environments
The pipeline was tested with the snakemake v7.32.4. The three conda environments on which the pipeline depends was deposited as yaml configure files in conda_environments directory. Users could create the corresponding environments using:
```
conda env create -f enrironment.yaml
```

## Modify the workflow configure file
All only thing the users need to do is modify the workflow configure file `config.yaml` to that located the same directory with snakemake workflow file. 
The `config.yaml` file include three types of parameters to set
1. Biosoftware path in your system
2. Parameters for some analysis steps
3. Input files path
   1. Users can get clear insight of the format for the most of required input data by scaning the example `config.yaml` file, such as fastq files of ATAC-seq data, fasta format of reference genome file
   2. For some custom input file, we provided example file in `annoFiles` directory
   
## Run workflow
After preparing all input data and correctely set values in `config.yaml` file, users can run the analysis workflow in a very easy way:
```
snakemake --use-conda -s MFVD.snakefile
```

# Expected output
1. **Allele specific open chromatin result**:`output_dir/ASoC_result/all.txt`, which contain the following columns:
	1. ***CHROM***: Chromosome
	2. ***POS***: Position of SNP
	3. ***REF***: Ref allele of the SNP
	4. ***ALT***: Alternative allele of the SNP
	5. ***heterozygosis_sample_count***: Number of samples in which the SNP was tested as heterozygote 
	6. ***REF_count***: Number of reads mapped to Ref allele across all heterozygous samples
	7. ***ALT_count***: Number of reads mapped to Alternative allele across all heterozygous samples
	8. ***DP_count***: Sum of REF_count and ALT_count
	9. ***P_value***: Pvalue of Binomial test
	10. ***FDR***: FDR calulated with p.adjust R package
2. **TFi-QTL result:** `output_dir/tf_regulation/tfieQTL_result_merged/tfieQTL.all.top_eVariant_for_each_TF_target_pair.withGeneName`, which contain the following main columns
	1. ***TF_id***: Ensemble id of TF; 
	2. ***TF_name***; Gene symbol of the TF
	3. ***target_id***; Ensemble id of eGene corresponding to `SNP`
	4. ***target_name***: Gene symbol of eGene corresponding to `SNP`
	5. ***SNP***: eSNP id
	6. ***interactP***: Pvalue of interaction term in regression model

## Workflow test
Users can download the ATAC-seq data (GSE108301) we used in our paper that deposite in bioRxiv 'https://www.biorxiv.org/content/10.1101/2024.06.11.598401v2.article-metrics' to test the workflow. It will cost several days to get the all results when using 16 ATAC-seq data. 
