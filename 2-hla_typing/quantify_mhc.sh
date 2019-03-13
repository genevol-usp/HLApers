#!/bin/bash

sample=$1
CPUS=$2

# Extract MHC and Unmapped reads
outPrefix=/scratch/genevol/users/vitor/mappings/${sample}_
bam=${outPrefix}Aligned.sortedByCoord.out.bam
fq1=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_1.fastq.gz
fq2=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_2.fastq.gz
readsalign=${outPrefix}readsAligned.txt
readsunmap=${outPrefix}readsUnmapped.txt
readids=${outPrefix}readids.txt
mhccoords=`cat ../1-make_indices/mhc_coords.txt`

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

fqmhc1=${outPrefix}mhc_1.fq
fqmhc2=${outPrefix}mhc_2.fq

seqtk subseq $fq1 $readids > $fq12 
seqtk subseq $fq2 $readids > $fq22 


# Remap to supplemented index
indexDIR=/home/vitor/hlapers/1-make_indices/indices/star
bammhcPrefix=${outPrefix}MHC_

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
fasta=/home/vitor/hlapers/1-make_indices/data/gencode/gencode.v25.MHC.IMGT.transcripts.fa
bammhc=${bammhcPrefix}Aligned.out.bam
out=${outPrefix}MHC_quants

if [ -d "$out" ]; then
    rm -r $out
fi

salmon quant -t $fasta -l IU -a $bammhc -o $out -p $CPUS

#Extract up to top 5 HLA alleles
Rscript write_top_alleles.R $out/quant.sf $out/top5_alleles.tsv

persindex=${outPrefix}persindex
persfasta=$persindex/hla.fa

if [ -d "$persindex" ]; then
    rm -r $persindex
fi
    
mkdir $persindex

Rscript write_genotyped_alleles.R $out/top5_alleles.tsv $persfasta

#Requantify expression of the top5
salmon index -t $persfasta -i $persindex/salmon --type quasi -k 31

outtop5=${outPrefix}top5_quants

mkdir $outtop5

salmon quant -i $persindex/salmon -l IU -1 $fqmhc1 -2 $fqmhc2 -o $outtop5\
    -p $CPUS --writeMappings > $outtop5/mappings.sam

Rscript genotype_top5.R $outtop5/quant.sf $outtop5

#Remove reads from the winner alleles
mappings=$outtop5/mappings.sam
winners=$outtop5/winners.txt
readsWin=${outPrefix}readsWin.txt
readsNoWin=${outPrefix}readsNoWin.txt

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

fqnoWin1=${outPrefix}noWin_1.fq
fqnoWin2=${outPrefix}noWin_2.fq

seqtk subseq $fqmhc1 $readsNoWin > $fqnoWin1
seqtk subseq $fqmhc2 $readsNoWin > $fqnoWin2

#Requantify to see if winner alleles explain all the expression or if 
# 2nd allele is still relevant
outNoWin=${outPrefix}NoWin_quants

mkdir -p $outNoWin

# HERE INDEX SHOULD HAVE ONLY 1ST AND 2ND ALLELES

salmon quant -i $persindex/salmon -l IU -1 $fqnoWin1 -2 $fqnoWin2\
    -o $outNoWin -p $CPUS

#Final gentotypes and personalized index
Rscript final_genotype.R $outtop5/quant.sf $outNoWin/quant.sf $outNoWin

Rscript write_genotyped_alleles.R $outNoWin/genotypes.tsv ${outPrefix}index.fa
