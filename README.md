# Instructions

## Install:

- R v3.4+
    - packages from CRAN:
	+ data.table 1.10.4+
	+ devtools v1.13.0+
	+ tidyverse v1.1.1+
	+ Biostrings v2.44.0
    - package from GitHub:
	+ hlaseqlib v0.0.0.9000+ (install with devtools::install\_github("vitoraguiar/hlaseqlib") inside R

- RSEM

- STAR v2.5.3a+

- Salmon v0.8.2+

- kallisto v0.43.1+

## Download data:

### IMGT

*scripts assume IMGT directory cloned in home directory*

Currently, the pipeline is tested with IMGT v3.29.0. To download this version:

git clone -b 3290 https://github.com/ANHIG/IMGTHLA.git 

Or, to download latest version:

git clone https://github.com/ANHIG/IMGTHLA.git


## The pipeline

### Building the index

The first step is to build an index composed of Gencode v25 transcripts, where
we replace the HLA isoforms with IMGT HLA allele sequences.

cd ./1-index_preparation

```
project
```
