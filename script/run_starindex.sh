#!/bin/bash

transcripts=$1
out=$2
threads=$3

mkdir -p $out

Rscript ./script/calc_starindex_params.R $transcripts $out

binbits=$(awk 'FNR == 1' $out/indexparams.txt)
saindex=$(awk 'FNR == 2' $out/indexparams.txt)

STAR --runThreadN $threads --runMode genomeGenerate \
    --genomeDir $out --genomeFastaFiles $transcripts \
    --genomeChrBinNbits $binbits --genomeSAindexNbases $saindex
