# Instructions

## Install:

- R v3.4+

- packages from CRAN:
    + data.table 1.10.4+
    + devtools v1.13.0+
    + tidyverse v1.1.1+
    + Biostrings v2.44.0

- package from GitHub:
    + hlaseqlib v0.0.0.9000+ 
    

- RSEM

- kallisto v0.43.1+

- STAR v2.5.3a+

- Salmon v0.8.2+


Notes:

- The hlaseqlib R package is under development. It can be installed with the
  command devtools::install\_github("vitoraguiar/hlaseqlib") inside R. As an
  alternative, if the user want to get the latest updates, we advise the user to
  clone the repository and, instead of loading the package with the standard
  "library()" function in R, use "devtools::load\_all(path\_to\_directory)".

- Users can choose between the faster kallisto and the STAR-Salmon pipeline.
  If STAR-Salmon is your choice, there is no need to install kallisto. 


## Download data:

### IMGT

*scripts assume IMGT directory cloned in home directory*

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

This script is not efficient, and if the users are allowed by their clusters to
use more CPUs when executing this script, we advise to modify the n\_cores
argument on line 14 to something like:

```
map2(locus, infer, hla_make_sequences, n_cores = 16)
```


Finally, we execute:

```
Rscript 4-make_index_fasta.R
```

This script creates 2 fasta files that will be used by the aligners to build the
indices. 

