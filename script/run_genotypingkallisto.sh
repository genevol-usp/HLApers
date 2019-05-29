#!/bin/bash

index=$1
transcripts=$2
fq1=$3
fq2=$4
out=$5
cpus=$6

kallisto quant -i $index -t $cpus -o $out $fq1 $fq2

Rscript genotype_kallisto.R $transcripts $out
