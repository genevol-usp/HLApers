#!/bin/bash

readarray -t samples < ../data/sample_ids.txt
quant=./quantifications_MHC
out=$quant/imgt_quants.tsv

awk 'FNR==1 {print "subject\t" $0}' $quant/${samples[0]}/quant.sf > $out

for id in "${samples[@]}"
do
    file=$quant/$id/quant.sf

    awk -v SUBJECT="$id" 'FNR>1 && $1 ~ /IMGT/ {print SUBJECT "\t" $0}' $file >> $out
    
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $out > $out.tmp && mv $out.tmp $out
