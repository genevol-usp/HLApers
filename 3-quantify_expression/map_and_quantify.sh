#!/bin/bash

STAR=~/STAR
salmon=~/Salmon-latest_linux_x86_64/bin/salmon

sample=$1
CPUS=6
indexDIR=./sample_indices/$sample
fasta=../1-make_indices/gencode.v25.PRI.transcripts.noIMGT.fa
sample_hla=../3-hla_typing/sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa

mkdir -p $indexDIR

cat $fasta $sample_hla > $sample_fa

$STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir $indexDIR\
    --genomeFastaFiles $sample_fa\
    --genomeChrBinNbits 11 --genomeSAindexNbases 13\
    --outFileNamePrefix ${indexDIR}_

fq1=../data/fastq/${sample}_1.fastq.gz
fq2=../data/fastq/${sample}_2.fastq.gz
outMap=./mappings
outPrefix=$outMap/${sample}_

$STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
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
outQuant=./quantifications
out=$outQuant/$sample

if [ -d "$out" ]; then
    rmdir $out
fi

$salmon quant -t $sample_fa -l IU -a $bam -o $out -p $CPUS --seqBias --gcBias

rm -r $indexDIR ${indexDIR}_Log.out $sample_fa $outPrefix* 
