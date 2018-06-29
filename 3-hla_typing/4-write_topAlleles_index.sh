#!/bin/bash

write_custom_index=../1-make_indices/write_genotyped_alleles.R

outdir=./sample_indices

if [ -d "$outdir" ]; then
    rm -r $outdir
fi
    
mkdir -p $outdir

Rscript $write_custom_index ./quantifications_MHC/imgt_quants_topAlleles.tsv $outdir
