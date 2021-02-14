#!/bin/bash

DIR=$( dirname "$0" )
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
bampers=${outPrefix}_Aligned.out.bam
out=${outPrefix}_quant

mkdir -p $index
mkdir -p ${outPrefix}_log

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

rm -r $bampers $index $sample_transcripts $samplehla ${outPrefix}_SJ* 
