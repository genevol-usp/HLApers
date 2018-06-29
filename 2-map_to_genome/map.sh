#!/bin/bash

STAR=~/STAR
samtools=~/samtools-1.3.1/samtools
seqtk=~/seqtk/seqtk

sample=$1

CPUS=8
indexDIR=../1-make_indices/genome_index
fq1=../data/fastq/${sample}_1.fastq.gz
fq2=../data/fastq/${sample}_2.fastq.gz
outMap=./mappings
outPrefix=${outMap}/${sample}_

$STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterMultimapNmax 20\
  --outSAMtype BAM SortedByCoordinate\
  --outSAMunmapped Within KeepPairs\
  --quantTranscriptomeBan Singleend\
  --outFileNamePrefix $outPrefix

bamGenome=${outPrefix}Aligned.sortedByCoord.out.bam
readsalign=${outPrefix}readsAligned.txt
readsunmap=${outPrefix}readsUnmapped.txt
readids=${outPrefix}readids.txt

$samtools index $bamGenome $bamGenome.bai

$samtools view $bamGenome "chr6:29722000-33144000" |\
    cut -f1 |\
    sort |\
    uniq > $readsalign

$samtools view -F 0x2 $bamGenome |\
    cut -f1 |\
    sort |\
    uniq > $readsunmap

cat $readsalign $readsunmap |\
    sort |\
    uniq > $readids

fq12=./mappings/mhc_fqs/${sample}_1.fq
fq22=./mappings/mhc_fqs/${sample}_2.fq

$seqtk subseq $fq1 $readids > $fq12 
$seqtk subseq $fq2 $readids > $fq22 

rm ${outPrefix}*
