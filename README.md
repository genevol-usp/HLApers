# Instructions

## Install:

- R v3.4+

- R packages from CRAN:
    + devtools v1.13.0+
    + tidyverse v1.1.1+
    
- R packages from Bioconductor:
    + Biostrings v2.44.0

- R packages from GitHub:
    + hlaseqlib v0.0.0.9000+ (https://github.com/genevol-usp/hlaseqlib) 
    
- RSEM (git clone https://github.com/deweylab/RSEM.git)

- STAR v2.5.3a+

- Salmon v0.8.2+


Notes:

- The hlaseqlib R package is under active development. It can be installed with
  the command devtools::install\_github("vitoraguiar/hlaseqlib") inside R. As an
  alternative, in order to get the latest updates, we advise the user to clone
  the repository and, instead of loading the package with the standard
  "library()" function in R, use "devtools::load\_all(path\_to\_directory)".

- Scripts assume that the executables or directories of the programs above are
  in the home directory. Users must adjust the path in each script, or add the
  programs to $PATH.

## Download data:

*scripts assume IMGT repository cloned in home directory*

```
git clone https://github.com/ANHIG/IMGTHLA.git
```

## The pipeline

The pipeline is under development. We provide a mix of scripts to be
submitted to a PBS qsub cluster, and scripts to be executed on the front end. 

The scripts are located on the directories where they are executed. Scripts are
named according to the order of execution in each directory. When script names
do not start with a number, the script is not meant to be executed by the user,
but it is called by other scripts.

This may change in a near feature.


### Building a transcriptome supplemented with HLA sequences

The first step is to build an index composed of Gencode v25 transcripts, where
we replace the HLA transcripts with IMGT HLA allele sequences.

```
cd ./1-make_indices
```

By executing 

```
./0-download_data.sh
```

we create the subdirecty "gencode", and download the annotations and Genome
sequence for the primary assembly of the reference genome (GRCh38).

By executing the scripts numbered 1--9, we create the fasta and index files
which will be used by the aligners.


### Extract HLA reads from a standard alignment to the Genome

Here we run a standard mapping to the reference genome, and extract the reads
mapping to the MHC, and also those which are unmapped.

If users already have a BAM file, it is possible to adapt the scripts to skip
read mapping and proceed to the read extraction.

Assuming that there is a directory ./data/fastq which contains the fastq files
for all individuals, we can execute the mapping script

```
cd 2-map_to_genome
qsub 1-map.pbs
```

### HLA typing

Next, we move to the directory 3-hla\_typing

```
cd 3-hla_typing
```

and execute the scripts in the order which they are numbered to perform the HLA
typing.

### Quantification

Finally, we estimate expression.

```
cd 4-quantify_expression
```

Here we also execute the scripts in order.

