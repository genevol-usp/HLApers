#!/bin/bash

hladb=$1
samplehla=$2
fq1=$3
fq2=$4
outPrefix=$5
cpus=$6

transcriptsNoHLA=$hladb/gencode_noHLA.fa
index=${outPrefix}_index
sample_transcripts=${outPrefix}_transcripts.fa
out=${outPrefix}_quant

cat $transcriptsNoHLA $samplehla > $sample_transcripts

kallisto index -i $index $sample_transcripts

kallisto quant -i $index -t $cpus -o $out --bias $fq1 $fq2

awk 'FNR == 1 {print $1"\t"$4"\t"$5}' $out/abundance.tsv > ${outPrefix}_hlaquant.tsv
awk '/IMGT/ {print $1"\t"$4"\t"$5}' $out/abundance.tsv >> ${outPrefix}_hlaquant.tsv

rm $sample_transcripts $index
