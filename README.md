# Instructions

## Install:

- R v3.4+

- R packages from CRAN:
    + data.table 1.10.4+
    + devtools v1.13.0+
    + tidyverse v1.1.1+
    
- R packages from Bioconductor:
    + Biostrings v2.44.0

- R packages from GitHub:
    + hlaseqlib v0.0.0.9000+ (https://github.com/genevol-usp/hlaseqlib) 
    
- RSEM (git clone https://github.com/deweylab/RSEM.git)

- kallisto v0.43.1+

- STAR v2.5.3a+

- Salmon v0.8.2+


Notes:

- The hlaseqlib R package is under active development. It can be installed with
  the command devtools::install\_github("vitoraguiar/hlaseqlib") inside R. As an
  alternative, in order to get the latest updates, we advise the user to clone
  the repository and, instead of loading the package with the standard
  "library()" function in R, use "devtools::load\_all(path\_to\_directory)".

- Users can choose between the faster kallisto and the STAR-Salmon pipeline.
  If STAR-Salmon is your choice, there is no need to install kallisto. 

- Scripts assume that the executables or directories of the programs above are
  in the home directory. Users must examine each script and adjust the path, or
  add the programs to the $PATH.

## Download data:

*scripts assume IMGT repository cloned in home directory*

Currently, the pipeline is tested with IMGT v3.29.0. To download this version:

```
git clone -b 3290 https://github.com/ANHIG/IMGTHLA.git
```

Or, to download latest version:

```
git clone https://github.com/ANHIG/IMGTHLA.git
```

## The pipeline

The pipeline is still under development. We provide a mix of scripts to be
submitted to a PBS qsub cluster and scripts to be executed outside the job
submission system (e.g., on the front end). 

Users outside the LEM-darwin cluster must adapt the pipeline if:

1. Their cluster uses a job submission system different than PBS qsub.

2. Their cluster does not allow scripts with potential high RAM usage to run on
   the front end.

The scripts are located on the directories where they are executed. This may
change in a near feature.

Scripts are named according to the order of execution in each directory. When 
script names do not start with a number, the script is not meant to be executed
by the user, but it is called by other scripts.


### Building a transcriptome supplemented with HLA sequences

The first step is to build an index composed of Gencode v25 transcripts, where
we replace the HLA isoforms with IMGT HLA allele sequences.

```
cd ./1-index_preparation
```

By executing 

```
./1-download_data.sh
```

we create the subdirecty "gencode", and download the annotations and Genome
sequence for the primary assembly of the reference genome (GRCh38).

Then we run RSEM to slice the reference genome according to the annotations
in order to generate a fasta file with transcript sequences. 

```
qsub 2-rsem_prepare_PRI.pbs
```

Next, we build a fasta with HLA sequences by processing the IMGT alignments.

```
Rscript 3-make_imgt_index.R
```

In this script, we hard code the names of HLA genes for which we want to infer 
exon sequences which are missing in IMGT.

line nº 5:
```
main_loci <- c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB")
```

The user may want to modify that list. Genes not listed here only contribute
with the alleles for which the complete sequence is available from IMGT.

This script is very inefficient (it takes 2h20 on our machine with a single
core), and if the users are allowed by their clusters to use more CPUs when
executing this script, we advise to modify the n\_cores argument on line 14 to
something like:

```
map2(locus, infer, hla_make_sequences, n_cores = 16)
```

Finally, we execute:

```
Rscript 4-make_index_fasta.R
```

This script creates 2 fasta files that will be used by the aligners to build the
indices. 


### Quantification

We will describe the STAR-Salmon pipeline, but the kallisto pipeline is similar
in structure.

```
cd 2-star_pipeline
```

The next script must be executed when all the jobs from the previous script are
finished. Or a pipeline.sh script can be made adding dependencies according to
your cluster (i.e., execute script 2 afterok jobs from script 1).


1. Build STAR index:

```
qsub 0-star_indexer.pbs
```

2. Run round 1 of quantification:

```
qsub 1-map_and_quantify_round1.pbs
```

This script assumes a text file in "hla\_rnaseq\_pipeline/data/sample\_ids.txt"
with sample IDs (one per line).

In this example script, we assume 100 samples. Modify the line 7 to adjust that.

The executable "./map\_and\_quantify.sh" is called. It assumes fastq files
located on "hla\_rnaseq\_pipeline/data/" and named as "SAMPLEID\_1.fastq.gz" and
"SAMPLEID\_2.fastq.gz".


3. When the quantification of all samples is finished, execute the following
   script:

```
./2-compile_quants_round1.sh
```

4. We then process the HLA quantifications and infer the HLA types:

```
Rscript 3-process_imgt_quants_round1.R
```

5. Write a custom index for each sample:

```
./4-write_custom_index.sh
```

6. Run round 2 of quantification:

```
qsub 5-map_and_quantify_round2.pbs
```

This script, and the executable "map\_and\_quantify\_round1.sh", must be edited
just like the version of these scripts for the round 1 to suit the number of
samples, fastq names etc.

7. After all quantification jobs are finished, compile the results:

```
./6-compile_quants_round2.sh
```

8. Process the HLA quantifications and infer HLA types:

```
Rscript 7-process_imgt_quants_round2.R
```

9. Finally, if you wish to write a BED file with the quantifications of genes in
the reference autosome chromosomes, execute:

```
Rscript 8-write_quantifications_bed.R
```

Genes located on scaffolds, such as HLA-DRB3 and HLA-DRB4 will no be present in
this BED file, unless the script is modified to include them.
