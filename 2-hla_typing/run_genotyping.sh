#!/bin/bash

sample=$1
CPUS=$2

outPrefix=/scratch/genevol/users/vitor/mappings/${sample}_
bam=${outPrefix}Aligned.sortedByCoord.out.bam
fq1=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_1.fastq.gz
fq2=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_2.fastq.gz
gencode=/home/vitor/hlapers/1-make_indices/data/gencode/gencode.v25.MHC.IMGT.transcripts.fa
readsalign=${outPrefix}readsAligned.txt
readsunmap=${outPrefix}readsUnmapped.txt
readids=${outPrefix}readids.txt
mhccoords=$(cat ../1-make_indices/mhc_coords.txt)
fqmhc1=${outPrefix}mhc_1.fq
fqmhc2=${outPrefix}mhc_2.fq
indexDIR=/home/vitor/hlapers/1-make_indices/indices/star
bammhcPrefix=${outPrefix}MHC_
bammhc=${bammhcPrefix}Aligned.out.bam
outmhc=${bammhcPrefix}quants
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

# Extract MHC and Unmapped reads
samtools index $bam $bam.bai

samtools view $bam $mhccoords |\
    cut -f1 |\
    sort |\
    uniq > $readsalign

samtools view -F 0x2 $bam |\
    cut -f1 |\
    sort |\
    uniq > $readsunmap

cat $readsalign $readsunmap |\
    sort |\
    uniq > $readids

seqtk subseq $fq1 $readids > $fqmhc1 
seqtk subseq $fq2 $readids > $fqmhc2 

# Remap to supplemented index
STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
    --readFilesIn $fqmhc1 $fqmhc2\
    --outFilterMismatchNmax 1\
    --outFilterMultimapScoreRange 0\
    --outFilterMultimapNmax 3000\
    --winAnchorMultimapNmax 6000\
    --alignEndsType EndToEnd\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix $bammhcPrefix

# Quantify MHC expression
if [ -d "$outmhc" ]; then
    rm -r $outmhc
fi

salmon quant -t $gencode -l IU -a $bammhc -o $outmhc -p $CPUS

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
    -p $CPUS --writeMappings > $outtop5/mappings.sam

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
    -o $outNoWin -p $CPUS

#Final gentotypes and personalized index
Rscript final_genotype.R $outtop5/quant.sf $outNoWin/quant.sf ${outPrefix}genotypes.tsv 
Rscript write_genotyped_alleles.R ${outPrefix}genotypes.tsv ${outPrefix}index.fa

rm -r $readids $readsalign $readsunmap $fqmhc1 $fqmhc2 $bammhc $outmhc \
    ${bammhcPrefix}SJ.out.tab $persindex $outtop5 $readsWin $readsNoWin \
    $fqnoWin1 $fqnoWin2 $outNoWin $bam $bam.bai
