#!/bin/bash

fileDir=$1
ic_threshold=$2
frequency_threshold=$3

for motifDatabase in HOCOMOCO jaspar
do
    inputFile="${fileDir}/variant.TF_footprint.overlap.${motifDatabase}.motif_disruption"
    awk -v v1="$ic_threshold" -v v2="$frequency_threshold" '{if(NR==1 || ($23>=v2 && $18>=v1)) print $0}' ${inputFile} > ${inputFile}.match_threshold
done