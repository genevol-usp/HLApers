#!/bin/bash

indexDIR=$1
sample=$2
fqDIR=$3
mhccoords=$4
gencode=$5
outDIR=$6
cpus=$7

outPrefix=$outDIR/${sample}_
bam=${outPrefix}Aligned.sortedByCoord.out.bam
fq1=$fqDIR/${sample}_1.fastq.gz
fq2=$fqDIR/${sample}_2.fastq.gz
readsalign=${outPrefix}readsAligned.txt
readsunmap=${outPrefix}readsUnmapped.txt
readids=${outPrefix}readids.txt
fqmhc1=${outPrefix}mhc_1.fq
fqmhc2=${outPrefix}mhc_2.fq
mhcPrefix=${outPrefix}MHC_
bammhc=${mhcPrefix}Aligned.out.bam
outmhc=${mhcPrefix}quants
persindex=${outPrefix}persindex
persfasta=$persindex/hla.fa
outtop5=${outPrefix}top5_quants
mappings=$outtop5/mappings.sam
winners=$outtop5/winners.txt
readsWin=${outPrefix}readsWin.txt
readsNoWin=${outPrefix}readsNoWin.txt
fqnoWin1=${outPrefix}noWin_1.fq
fqnoWin2=${outPrefix}noWin_2.fq
outNoWin=${outPrefix}NoWin_quants

## Extract MHC and Unmapped reads
#samtools index $bam $bam.bai
#
#samtools view $bam $mhccoords |\
#    cut -f1 |\
#    sort |\
#    uniq > $readsalign
#
#samtools view -F 0x2 $bam |\
#    cut -f1 |\
#    sort |\
#    uniq > $readsunmap
#
#cat $readsalign $readsunmap |\
#    sort |\
#    uniq > $readids
#
#seqtk subseq $fq1 $readids > $fqmhc1 
#seqtk subseq $fq2 $readids > $fqmhc2 
#
## Remap to supplemented index
STAR --runMode alignReads --runThreadN $cpus --genomeDir $indexDIR\
    --readFilesIn $fqmhc1 $fqmhc2\
    --outFilterMismatchNmax 1\
    --outFilterMultimapScoreRange 0\
    --outFilterMultimapNmax 3000\
    --winAnchorMultimapNmax 6000\
    --alignEndsType EndToEnd\
    --alignTranscriptsPerReadNmax 100000\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix $mhcPrefix

# Quantify MHC expression
if [ -d "$outmhc" ]; then
    rm -r $outmhc
fi

salmon quant -t $gencode -l IU -a $bammhc -o $outmhc -p $cpus

#Extract up to top 5 HLA alleles
if [ -d "$persindex" ]; then
    rm -r $persindex
fi
    
mkdir $persindex
mkdir $outtop5

Rscript write_top_alleles.R $outmhc/quant.sf $outtop5/top5_alleles.tsv
Rscript write_genotyped_alleles.R $outtop5/top5_alleles.tsv $persfasta

#Requantify expression of the top5
salmon index -t $persfasta -i $persindex/salmon --type quasi -k 31

salmon quant -i $persindex/salmon -l IU -1 $fqmhc1 -2 $fqmhc2 -o $outtop5\
    -p $cpus --writeMappings > $outtop5/mappings.sam

Rscript genotype_top5.R $outtop5/quant.sf $outtop5

#Remove reads from the winner alleles
grep -v "^@" $mappings |\
    grep -F -f $winners - |\
    cut -f1 |\
    sort |\
    uniq > $readsWin

grep -v "^@" $mappings |\
    cut -f1 |\
    awk 'FNR==NR {hash[$0]; next} !($0 in hash)' $readsWin - |\
    sort |\
    uniq > $readsNoWin

seqtk subseq $fqmhc1 $readsNoWin > $fqnoWin1
seqtk subseq $fqmhc2 $readsNoWin > $fqnoWin2

#Requantify to see if winner alleles explain all the expression or if 
# there is other relevant allele
mkdir -p $outNoWin

salmon quant -i $persindex/salmon -l IU -1 $fqnoWin1 -2 $fqnoWin2\
    -o $outNoWin -p $cpus

#Final gentotypes and personalized index
Rscript final_genotype.R $outtop5/quant.sf $outNoWin/quant.sf ${outPrefix}genotypes.tsv 
Rscript write_genotyped_alleles.R ${outPrefix}genotypes.tsv ${outPrefix}index.fa

mv ${mhcPrefix}Log* $outDIR/log/ 
mv $outmhc/logs/salmon_quant.log $outDIR/log/${sample}_salmon_quant_mhc.log

#rm -r $readids $readsalign $readsunmap $fqmhc1 $fqmhc2 $bammhc $outmhc \
#    ${mhcPrefix}SJ.out.tab $persindex $outtop5 $readsWin $readsNoWin \
#    $fqnoWin1 $fqnoWin2 $outNoWin 
#
#
rm -r $bammhc $outmhc\
    ${mhcPrefix}SJ.out.tab $persindex $outtop5 $readsWin $readsNoWin \
    $fqnoWin1 $fqnoWin2 $outNoWin 
