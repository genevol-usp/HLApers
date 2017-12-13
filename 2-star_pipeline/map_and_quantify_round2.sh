#!/bin/bash

STAR=~/STAR
salmon=~/Salmon-latest_linux_x86_64/bin/salmon

sample=$1
indexDIR=./sample_indices/$sample
sample_hla=./sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa
fasta=../1-index_preparation/gencode.v25.PRI.transcripts.noIMGT.fa

mkdir -p $indexDIR

cat $fasta $sample_hla > $sample_fa

$STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $indexDIR\
    --genomeFastaFiles $sample_fa\
    --genomeChrBinNbits 11 --genomeSAindexNbases 13\
    --outFileNamePrefix ${indexDIR}_

fq1=../data/fastq/${sample}_1.fastq.gz
fq2=../data/fastq/${sample}_2.fastq.gz
outMap=./mappings_2
outQuant=./quantifications_2
outPrefix=$outMap/${sample}_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
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
    --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.out.bam
out=$outQuant/$sample

if [ -d "$out" ]; then
    rm -r $out
fi

$salmon quant -t $sample_fa -l IU -a $bam -o $out -p 6 --seqBias --gcBias

rm -r $indexDIR ${indexDIR}_Log.out $sample_fa $sample_hla $outPrefix* 
