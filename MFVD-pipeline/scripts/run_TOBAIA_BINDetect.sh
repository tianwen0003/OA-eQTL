#!/bin/bash

tobiasPath=$1
threads=$2
bw_signal=$3
ref=$4
alter_ref=$5
peaks=$6
jasparMotif=$7
hocomodoMotif=$8
outputDir=$9

#----jaspar using reference sequence----#
${tobiasPath} BINDetect --signal ${bw_signal} \
                     --genome ${ref} \
                     --peaks ${peaks} \
                     --motifs ${jasparMotif} \
                     --outdir ${outputDir}/03_BINDetect_jaspar_originalRef \
                     --prefix all_samples \
                     --cores ${threads}

#----jaspar using alternative allele replaced reference sequence----#
${tobiasPath} BINDetect --signal ${bw_signal} \
                     --genome ${alter_ref} \
                     --peaks ${peaks} \
                     --motifs ${jasparMotif} \
                     --outdir ${outputDir}/03_BINDetect_jaspar_snpsIntroducedRef \
                     --prefix all_samples \
                     --cores ${threads}

#----HOCOMOCO using reference sequence----#
${tobiasPath} BINDetect --signal ${bw_signal} \
                     --genome ${ref} \
                     --peaks ${peaks} \
                     --motifs ${hocomodoMotif} \
                     --outdir ${outputDir}/03_BINDetect_HOCOMOCO_originalRef \
                     --prefix all_samples \
                     --cores ${threads}

#----jaspar using alternative allele replaced reference sequence----#
${tobiasPath} BINDetect --signal ${bw_signal} \
                     --genome ${alter_ref} \
                     --peaks ${peaks} \
                     --motifs ${hocomodoMotif} \
                     --outdir ${outputDir}/03_BINDetect_HOCOMOCO_snpsIntroducedRef \
                     --prefix all_samples \
                     --cores ${threads}
