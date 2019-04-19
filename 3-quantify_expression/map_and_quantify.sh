#!/bin/bash

sample=$1
CPUS=$2

genocodeNoIMGT=../1-make_indices/data/gencode/gencode.v25.PRI.transcripts.noIMGT.fa
prefix=/scratch/genevol/users/vitor/mappings/${sample} 
indexDIR=${prefix}_index
sample_hla=${prefix}_index.fa
sample_gencode=${prefix}_gencode.fa
fq1=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_1.fastq.gz
fq2=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_2.fastq.gz

#mkdir -p $indexDIR
#
#cat $genocodeNoIMGT $sample_hla > $sample_gencode
#
#STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir $indexDIR\
#    --genomeFastaFiles $sample_gencode\
#    --genomeChrBinNbits 11 --genomeSAindexNbases 13\
#    --outFileNamePrefix ${indexDIR}_
#

STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
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
    --outFileNamePrefix ${prefix}_

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
