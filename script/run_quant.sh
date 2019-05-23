#!/bin/bash

bam=$1
hladb=$2
samplehla=$3
outPrefix=$4
cpus=$5

gencodeNoHLA=$hladb/gencode_noHLA.fa
index=${outPrefix}_index
sample_gencode=${outPrefix}_gencode.fa
fq1=${outPrefix}_1.fq.gz
fq2=${outPrefix}_2.fq.gz
bampers=${outPrefix}_Aligned.out.bam
out=${outPrefix}_quant

samtools sort -n $bam |\
    samtools fastq -N -1 $fq1 -2 $fq2 -

mkdir -p $index

cat $gencodeNoHLA $samplehla > $sample_gencode

STAR --runThreadN $cpus --runMode genomeGenerate --genomeDir $index\
    --genomeFastaFiles $sample_gencode\
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
    --outSAMunmapped Within KeepPairs\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix ${outPrefix}_

salmon quant -t $sample_gencode -l A -a $bampers -o $out -p $cpus \
    --seqBias --gcBias --posBias

mv ${outPrefix}_*Log.* ${outPrefix}_logs/

rm -r $bampers $fq1 $fq2 $index $sample_gencode ${outPrefix}_SJ* 
