#!/bin/bash

# create dir gencode
mkdir -p ./gencode

# download genome and annotations
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.primary_assembly.genome.fa.gz -P ./gencode
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.primary_assembly.annotation.gtf.gz -P ./gencode

# download IMGT database
git clone https://github.com/ANHIG/IMGTHLA.git

# create transcript fasta from genome and annotations with RSEM
genome=./gencode/GRCh38.primary_assembly.genome.fa
gtf=./gencode/gencode.v25.primary_assembly.annotation.gtf
out=./gencode/gencode.v25.PRI

zcat $genome.gz > $genome
zcat $gtf.gz > $gtf

rsem-prepare-reference --gtf $gtf $genome $out

rm $out.grp $out.ti $out.chrlist $out.seq $out.idx.fa $out.n2g.idx.fa $genome $gtf
