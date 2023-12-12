#!/bin/bash

bedtools=$1
variants_file=$2
input_dir=$3
output_dir=$4

for tfDatabase in HOCOMOCO jaspar
do  
    sed '1d' ${variants_file} \
    | ${bedtools} intersect -wo -a - -b ${input_dir}/${tfDatabase}.all_motif_footprint.bed \
    | cut -f 1-17 \
    | sed '1i eVariant_chromosome\tvariant_position_start\tvariant_position_end\tvariant_id\teQTL_gene_id\teQTL_gene_symbol\teffect_allele\tbaseLine_allele\tmaf\tslope\tmotif_chromosome\tmotif_start\tmotif_end\tmotif_id\tTF_gene_symbol\tTF_gene_id\tmotif_strand' \
    > ${output_dir}/variant.TF_footprint.overlap.${tfDatabase}
done