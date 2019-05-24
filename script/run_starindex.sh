#!/bin/bash

transcripts=$1
out=$2
threads=$3

mkdir -p $out/STARMHC

STAR --runThreadN $threads --runMode genomeGenerate \
    --genomeDir $out/STARMHC --genomeFastaFiles $transcripts \
    --genomeChrBinNbits 10 --genomeSAindexNbases 11

