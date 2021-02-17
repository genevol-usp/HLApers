#!/bin/bash

DIR=$( dirname "$0" )
index=$1
fq1=$2
fq2=$3
out=$4
cpus=$5

kallisto quant -i $index -t $cpus -o $out $fq1 $fq2

Rscript $DIR/genotype_kallisto.R $out
