#!/bin/bash

bam=$1
index=$2
gencode=$3
mhccoords=$4
outPrefix=$5
cpus=$6

mhc=$(cat $mhccoords)

#bam to fastq
mhcfq1=${outPrefix}_mhc_1.fq 
mhcfq2=${outPrefix}_mhc_2.fq 
unmapfq1=${outPrefix}_unmap_1.fq
unmapfq2=${outPrefix}_unmap_2.fq
samplefq1=${outPrefix}_1.fq
samplefq2=${outPrefix}_2.fq

# remapping
bammhc=${outPrefix}_MHC_Aligned.out.bam
outmhc=${outPrefix}_MHC_quants
persindex=${outPrefix}_persindex
outtop5=${outPrefix}_top5_quants
readsWin=${outPrefix}_readsWin.txt
readsNoWin=${outPrefix}_readsNoWin.txt
fqnoWin1=${outPrefix}_noWin_1.fq
fqnoWin2=${outPrefix}_noWin_2.fq
outNoWin=${outPrefix}_NoWin_quants

# Extract MHC reads

echo "Extracting MHC and unmapped reads from BAM..."

if [ ! -f "$bam".bai ]; then
    samtools index $bam $bam.bai
fi

samtools view $bam $mhc -b |\
    samtools sort -n - |\
    samtools fastq -1 $mhcfq1 -2 $mhcfq2 -

samtools view -F 0x2 $bam |\
    samtools sort -n - |\
    samtools fastq -1 $unmapfq1 -2 $unmapfq2 -

cat $mhcfq1 $unmapfq1 > $samplefq1.tmp
cat $mhcfq2 $unmapfq2 > $samplefq2.tmp

comm -12 <(sed -n '1~4p' $samplefq1.tmp | sort) <(sed -n '1~4p' $samplefq2.tmp | sort) |\
    sed 's|^@||' |\
    sort -V |\
    uniq > ${outPrefix}_reads

seqtk subseq $samplefq1.tmp ${outPrefix}_reads > $samplefq1
seqtk subseq $samplefq2.tmp ${outPrefix}_reads > $samplefq2

rm $mhcfq1 $mhcfq2 $unmapfq1 $unmapfq2 $samplefq1.tmp $samplefq2.tmp ${outPrefix}_reads

# Remap to supplemented index

echo "Remapping extracted reads to personalized MHC index..."

STAR --runMode alignReads --runThreadN $cpus --genomeDir $index\
    --readFilesIn $samplefq1 $samplefq2\
    --outFilterMismatchNmax 1\
    --outFilterMultimapScoreRange 0\
    --outFilterMultimapNmax 3000\
    --winAnchorMultimapNmax 6000\
    --alignEndsType EndToEnd\
    --alignTranscriptsPerReadNmax 100000\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix ${outPrefix}_MHC_  

## Quantify MHC expression
#salmon quant -t $gencode -l IU -a $bammhc -o $outmhc -p $cpus
#
#Extract up to top 5 HLA alleles
#mkdir -p $persindex
#
#Rscript ./script/write_top5_fasta.R $outmhc/quant.sf $persindex/hla.fa
#
##Requantify expression of the top5
#mkdir -p $outtop5
#
#salmon index -t $persindex/hla.fa -i $persindex/salmon --type quasi -k 31
#
#salmon quant -i $persindex/salmon -l IU -1 $samplefq1 -2 $samplefq2 -o $outtop5\
#    -p $cpus --writeMappings > $outtop5/mappings.sam
#
#Rscript ./script/write_winners.R $outtop5/quant.sf $outtop5/winners.txt
#
##Remove reads from the winner alleles
#grep -v "^@" $outtop5/mappings.sam |\
#    grep -F -f $outtop5/winners.txt - |\
#    cut -f1 |\
#    sort |\
#    uniq > $readsWin
#
#grep -v "^@" $outtop5/mappings.sam|\
#    cut -f1 |\
#    awk 'FNR==NR {hash[$0]; next} !($0 in hash)' $readsWin - |\
#    sort |\
#    uniq > $readsNoWin
#
#seqtk subseq $samplefq1 $readsNoWin > $fqnoWin1
#seqtk subseq $samplefq2 $readsNoWin > $fqnoWin2
#
#Requantify to see if winner alleles explain all the expression or if 
# there is other relevant allele
#
#salmon quant -i $persindex/salmon -l IU -1 $fqnoWin1 -2 $fqnoWin2\
#    -o $outNoWin -p $cpus
#
##Final gentotypes and personalized index
#Rscript ./script/write_final_genotypes.R $outtop5/quant.sf $outNoWin/quant.sf ${outPrefix} 
#
#mkdir -p ${outPrefix}_logs
#
#mv ${outPrefix}_MHC_Log* ${outPrefix}_logs/
#mv ${outPrefix}_MHC_quants/logs/salmon_quant.log ${outPrefix}_logs/ 
#
#rm -r $samplefq1 $samplefq2 ${outPrefix}_MHC* $persindex $outtop5 $readsWin\
#    $readsNoWin $fqnoWin1 $fqnoWin2 $outNoWin 
#
