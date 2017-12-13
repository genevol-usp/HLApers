#!/bin/bash

quantDir=./quantifications_1
OUTimgt=$quantDir/imgt_quants.tsv

readarray -t samples < ../data/sample_ids.txt

awk 'FNR == 1 {print "subject\t" $0}' $quantDir/${samples[0]}/quant.sf > $OUTimgt.tmp

for id in "${samples[@]}" 
do
    file=$quantDir/$id/quant.sf

    awk -v SUBJECT="$id" 'FNR > 1 && $1 ~ /IMGT/ {print SUBJECT "\t" $0}' $file >> $OUTimgt.tmp
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $OUTimgt.tmp > $OUTimgt && rm $OUTimgt.tmp
