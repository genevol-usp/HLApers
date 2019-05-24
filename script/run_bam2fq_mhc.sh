#!/bin/bash

bam=$1
mhccoords=$2
outPrefix=$3

mhc=$(cat $mhccoords)

#bam to fastq
mapfq1=${outPrefix}_map_1.fq 
mapfq2=${outPrefix}_map_2.fq 
unmapfq1=${outPrefix}_unmap_1.fq
unmapfq2=${outPrefix}_unmap_2.fq
mhcfq1=${outPrefix}_mhc_1.fq
mhcfq2=${outPrefix}_mhc_2.fq

echo "Extracting MHC and unmapped reads from BAM..."

if [[ ! -f "$bam".bai ]]; then
    samtools index $bam $bam.bai
fi

samtools view $bam $mhc -b |\
    samtools sort -n - |\
    samtools fastq -1 $mapfq1 -2 $mapfq2 -

samtools view -F 0x2 $bam |\
    samtools sort -n - |\
    samtools fastq -1 $unmapfq1 -2 $unmapfq2 -

cat $mapfq1 $unmapfq1 > $mhcfq1.tmp
cat $mapfq2 $unmapfq2 > $mhcfq2.tmp

comm -12 <(sed -n '1~4p' $mhcfq1.tmp | sort) <(sed -n '1~4p' $mhcfq2.tmp | sort) |\
    sed 's|^@||' |\
    sort -V |\
    uniq > ${outPrefix}_reads

echo "Writing fastq files..."

seqtk subseq $mhcfq1.tmp ${outPrefix}_reads > $mhcfq1
seqtk subseq $mhcfq2.tmp ${outPrefix}_reads > $mhcfq2

rm $mapfq1 $mapfq2 $unmapfq1 $unmapfq2 $mhcfq1.tmp $mhcfq2.tmp ${outPrefix}_reads

echo "Done!"
