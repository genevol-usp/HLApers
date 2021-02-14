#!/bin/bash

DIR=$( dirname "$0" )
hladb=$1
genos=$2
fq1=$3
fq2=$4
outPrefix=$5
cpus=$6

out=${outPrefix}_quant
index=${outPrefix}_index
transcripts=$hladb/transcripts_MHC_HLAsupp.fa
transcriptsNoHLA=$hladb/transcripts_noHLA.fa
samplehla=${outPrefix}_index.fa
sample_transcripts=${outPrefix}_transcripts.fa

mkdir -p ${outPrefix}_log

Rscript $DIR/write_pers_index.R $transcripts $genos $outPrefix 

cat $transcriptsNoHLA $samplehla > $sample_transcripts

salmon index -t $sample_transcripts -i $index

salmon quant -i $index -l A -1 $fq1 -2 $fq2 -o $out -p $cpus \
    --seqBias --gcBias --posBias

awk 'FNR == 1 {print $1"\t"$4"\t"$5}' $out/quant.sf > ${outPrefix}_hlaquant.tsv
awk '/IMGT/ {print $1"\t"$4"\t"$5}' $out/quant.sf >> ${outPrefix}_hlaquant.tsv

rm -r $index $samplehla $sample_transcripts 
