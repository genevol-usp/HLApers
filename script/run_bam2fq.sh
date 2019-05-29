#!/bin/bash

bam=$1
outPrefix=$2

fq1=${outPrefix}_1.fq.gz
fq2=${outPrefix}_2.fq.gz

if [[ ! -f "$bam".bai ]]; then
    samtools index $bam $bam.bai
fi

samtools sort -n $bam |\
    samtools fastq -N -1 $fq1 -2 $fq2 -0 /dev/null -

