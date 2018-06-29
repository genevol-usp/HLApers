#!/bin/bash

zcat ./gencode/gencode.v25.primary_assembly.annotation.gtf.gz |\
    grep -F -f ./imgt_pri_ids.txt - > ./imgt.annotations.gtf

