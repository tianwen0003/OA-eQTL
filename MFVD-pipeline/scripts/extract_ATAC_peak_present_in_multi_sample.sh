#!/bin/bash

bedops=$1
bedtools=$2
sample_count=$3
input_dir=$4
output_dir=$5
mkdir -p ${output_dir}

for peakFile in ${@:6:$#};do
    idx=`basename ${peakFile}`
    cut -f 1-3 ${input_dir}/${idx} | awk -vidx=${idx} '{ print $0"\t"idx; }' > ${output_dir}/${idx}.id.bed
done

bedops --everything ${output_dir}/*id.bed | bedmap --echo --echo-map-id-uniq --delim '\t' - > ${output_dir}/ATAC_peak_all.bed
awk -vthreshold=${sample_count} '(split($5,ids,";") >= threshold)' ${output_dir}/ATAC_peak_all.bed > ${output_dir}/ATAC_peak_in_multiple_samples.bed
bedtools merge -i ${output_dir}/ATAC_peak_in_multiple_samples.bed > ${output_dir}/ATAC_peak_in_multiple_samples.merged.bed


# rm ${output_dir}/*id.bed
# rm ${output_dir}/ATAC_peak_all.bed
# rm ${output_dir}/ATAC_peak_in_multiple_samples.bed