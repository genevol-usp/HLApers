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
bampers=${outPrefix}_Aligned.out.bam
out=${outPrefix}_quant

mkdir -p $index
mkdir -p ${outPrefix}_log/

cat $transcriptsNoHLA $samplehla > $sample_transcripts

STAR --runThreadN $cpus --runMode genomeGenerate --genomeDir $index\
    --genomeFastaFiles $sample_transcripts\
    --genomeChrBinNbits 11 --genomeSAindexNbases 13\
    --outFileNamePrefix ${index}_

STAR --runMode alignReads --runThreadN $cpus --genomeDir $index\
    --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
    --outFilterMismatchNmax 999\
    --outFilterMismatchNoverReadLmax 0.04\
    --outFilterMultimapScoreRange 1\
    --outFilterMultimapNmax 150\
    --winAnchorMultimapNmax 300\
    --alignIntronMax 0\
    --alignEndsType Local\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix ${outPrefix}_

salmon quant -t $sample_transcripts -l A -a $bampers -o $out -p $cpus \
    --seqBias --gcBias --posBias

awk 'FNR == 1 {print $1"\t"$4"\t"$5}' $out/quant.sf > ${outPrefix}_hlaquant.tsv
awk '/IMGT/ {print $1"\t"$4"\t"$5}' $out/quant.sf >> ${outPrefix}_hlaquant.tsv

mv ${outPrefix}_*Log.* ${outPrefix}_log/

rm -r $bampers $index $sample_transcripts ${outPrefix}_SJ* 
