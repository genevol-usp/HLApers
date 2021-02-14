#!/bin/bash

DIR=$( dirname "$0")
hladb=$1
genos=$2
fq1=$3
fq2=$4
outPrefix=$5
cpus=$6

transcripts=$hladb/transcripts_MHC_HLAsupp.fa

Rscript $DIR/write_pers_index.R $transcripts $genos $outPrefix

transcriptsNoHLA=$hladb/transcripts_noHLA.fa
samplehla=${outPrefix}_index.fa
sample_transcripts=${outPrefix}_transcripts.fa

cat $transcriptsNoHLA $samplehla > $sample_transcripts

index=${outPrefix}_index
out=${outPrefix}_quant

kallisto index -i $index $sample_transcripts

kallisto quant -i $index -t $cpus -o $out --bias $fq1 $fq2

awk 'FNR == 1 {print $1"\t"$4"\t"$5}' $out/abundance.tsv > ${outPrefix}_hlaquant.tsv
awk '/IMGT/ {print $1"\t"$4"\t"$5}' $out/abundance.tsv >> ${outPrefix}_hlaquant.tsv

rm $samplehla $sample_transcripts $index
