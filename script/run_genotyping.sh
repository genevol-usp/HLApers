#!/bin/bash

DIR=$( dirname "$0" )
index=$1
transcripts=$2
fq1=$3
fq2=$4
outPrefix=$5
cpus=$6

bammhc=${outPrefix}_MHC_Aligned.out.bam
bammhcsort=${outPrefix}_MHC_Aligned_sorted.out.bam
outmhc=${outPrefix}_MHC_quants
persindex=${outPrefix}_persindex
outtop5=${outPrefix}_top5_quants
readsWin=${outPrefix}_readsWin.txt
readsNoWin=${outPrefix}_readsNoWin.txt
fqnoWin1=${outPrefix}_noWin_1.fq
fqnoWin2=${outPrefix}_noWin_2.fq
outNoWin=${outPrefix}_NoWin_quants
hlabed=$( dirname $transcripts )/hla.bed

if file $fq1 | grep -q gzip ; then
    readcommand=zcat
else
    readcommand="-"
fi

# Remap to supplemented index
echo "Remapping extracted reads to personalized MHC index..."

STAR --runMode alignReads --runThreadN $cpus --genomeDir $index\
    --readFilesIn $fq1 $fq2 --readFilesCommand $readcommand\
    --outFilterMismatchNmax 1\
    --outFilterMultimapScoreRange 0\
    --outFilterMultimapNmax 3000\
    --winAnchorMultimapNmax 6000\
    --alignEndsType EndToEnd\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix ${outPrefix}_MHC_  


samtools sort -@ $cpus -o $bammhcsort $bammhc
samtools index -@ $cpus $bammhcsort
samtools depth -a -b $hlabed -o ${outPrefix}.depth $bammhcsort

# Quantify MHC expression
echo "Genotyping HLA..." 

salmon quant -t $transcripts -l A -a $bammhc -o $outmhc -p $cpus

#Extract up to top 5 HLA alleles
mkdir -p $persindex

Rscript $DIR/write_top5_fasta.R $outmhc/quant.sf $transcripts ${outPrefix}.depth $persindex/hla.fa

#Requantify expression of the top5
mkdir -p $outtop5

salmon index -t $persindex/hla.fa -i $persindex/salmon

salmon quant -i $persindex/salmon -l A -1 $fq1 -2 $fq2 -o $outtop5\
    -p $cpus --writeMappings=$outtop5/mappings.sam

Rscript $DIR/write_winners.R $outtop5/quant.sf $outtop5/winners.txt

#Remove reads from the winner alleles
samtools view $outtop5/mappings.sam |\
    grep -F -f $outtop5/winners.txt - |\
    cut -f1 |\
    sort |\
    uniq > $readsWin

samtools view $outtop5/mappings.sam |\
    cut -f1 |\
    awk 'FNR==NR {hash[$0]; next} !($0 in hash)' $readsWin - |\
    sort |\
    uniq > $readsNoWin

if [[ $(head -n1 $fq1) =~ /1$ ]] && [[ $(head -n1 $fq2) =~ /2$ ]]; then

    awk '{ print $0 "/1" }' $readsNoWin > ${readsNoWin}1
    awk '{ print $0 "/2" }' $readsNoWin > ${readsNoWin}2

    seqtk subseq $fq1 ${readsNoWin}1 > $fqnoWin1
    seqtk subseq $fq2 ${readsNoWin}2 > $fqnoWin2

else

    seqtk subseq $fq1 $readsNoWin > $fqnoWin1
    seqtk subseq $fq2 $readsNoWin > $fqnoWin2

fi

#Requantify to see if winner alleles explain all the expression or if 
# there is other relevant allele

if [[ -s "$fqnoWin1" ]] && [[ -s "$fqnoWin2" ]]; then 
    salmon quant -i $persindex/salmon -l A -1 $fqnoWin1 -2 $fqnoWin2\
	-o $outNoWin -p $cpus
fi

#Final gentotypes
Rscript $DIR/write_final_genotypes.R $transcripts $outtop5/quant.sf $outNoWin/quant.sf $outPrefix 

mkdir -p ${outPrefix}_log

mv ${outPrefix}_MHC_Log* ${outPrefix}_log/
mv ${outPrefix}_MHC_quants/logs/salmon_quant.log ${outPrefix}_log/ 

#rm -r ${outPrefix}_MHC* $persindex $outtop5 $readsWin\
#    ${readsNoWin}* $fqnoWin1 $fqnoWin2 $outNoWin 
#
echo "Done!" 
