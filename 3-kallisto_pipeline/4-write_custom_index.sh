#!/bin/bash

write_custom_index=~/hla_rnaseq_pipeline/1-index_preparation/write_genotyped_alleles.R

outdir=./sample_indices

mkdir -p $outdir

Rscript $write_custom_index ./quantifications_1/processed_imgt_quants.tsv $outdir
