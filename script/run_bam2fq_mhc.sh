#!/bin/bash

bam=$1
mhccoords=$2
outPrefix=$3

mhc=$(cat $mhccoords)

mapbam=${outPrefix}_map.bam 
unmapbam=${outPrefix}_unmap.bam
tmpbam=${outPrefix}_tmp.bam
tmpfq1=${outPrefix}_tmp_1.fq
tmpfq2=${outPrefix}_tmp_2.fq
reads1tmp=${outPrefix}_reads1.tmp
reads2tmp=${outPrefix}_reads2.tmp
reads=${outPrefix}_reads
reads1=${outPrefix}_reads1
reads2=${outPrefix}_reads2
finalfq1=${outPrefix}_mhc_unmap_1.fq
finalfq2=${outPrefix}_mhc_unmap_2.fq

echo "Extracting MHC and unmapped reads from BAM..."

if [[ ! -f "$bam".bai ]]; then
    samtools index $bam $bam.bai
fi

samtools view $bam $mhc -b -o $mapbam
samtools view -F 0x2 $bam -b -o $unmapbam
samtools merge $tmpbam $mapbam $unmapbam 

samtools sort -n $tmpbam | samtools fastq -1 $tmpfq1 -2 $tmpfq2 -0 /dev/null -

sed -n '1~4p' $tmpfq1 | sed 's|^@||' | sed 's|/1$||'| sort > $reads1tmp
sed -n '1~4p' $tmpfq2 | sed 's|^@||' | sed 's|/2$||'| sort > $reads2tmp

comm -12 $reads1tmp $reads2tmp | sort -V | uniq > $reads

if [[ $(head -n1 $tmpfq1) =~ /1$ ]] && [[ $(head -n1 $tmpfq2) =~ /2$ ]]; then

    awk '{ print $0 "/1" }' $reads > $reads1
    awk '{ print $0 "/2" }' $reads > $reads2

else

    cp $reads $reads1
    cp $reads $reads2

fi

echo "Writing fastq files..."

seqtk subseq $tmpfq1 $reads1 > $finalfq1
seqtk subseq $tmpfq2 $reads2 > $finalfq2

rm $mapbam $unmapbam $tmpbam $tmpfq1 $tmpfq2 $reads1tmp $reads2tmp $reads $reads1 $reads2

echo "Done!"
