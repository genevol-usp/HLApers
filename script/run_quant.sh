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

samtools sort -n $bam |\
    samtools fastq -1 $fq1 -2 $fq2 -


#mkdir -p $indexDIR
#
#cat $gencodeNoIMGT $sample_hla > $sample_gencode
#
#STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir $indexDIR\
#    --genomeFastaFiles $sample_gencode\
#    --genomeChrBinNbits 11 --genomeSAindexNbases 13\
#    --outFileNamePrefix ${indexDIR}_
#
#
#STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
#    --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
#    --outFilterMismatchNmax 999\
#    --outFilterMismatchNoverReadLmax 0.04\
#    --outFilterMultimapScoreRange 1\
#    --outFilterMultimapNmax 150\
#    --winAnchorMultimapNmax 300\
#    --alignIntronMax 0\
#    --alignEndsType Local\
#    --outSAMunmapped Within KeepPairs\
#    --outSAMprimaryFlag AllBestScore\
#    --outSAMtype BAM Unsorted\
#    --outFileNamePrefix ${prefix}_
#
#bam=${outPrefix}Aligned.out.bam
#outQuant=./quantifications
#out=$outQuant/$sample
#
#if [ -d "$out" ]; then
#    rmdir $out
#fi
#
#$salmon quant -t $sample_fa -l IU -a $bam -o $out -p $CPUS --seqBias --gcBias
#
#rm -r $indexDIR ${indexDIR}_Log.out $sample_fa $outPrefix* 
