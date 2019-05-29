#!/bin/bash

transcripts=$1
out=$2

kallisto index -i $out $transcripts
