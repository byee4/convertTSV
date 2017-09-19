#!/usr/bin/env bash

python convertTSV.py \
--input ../data/filtered_gene_bc_matrices_h5.h5 \
--output ../data/outputs/matrix.tsv \
--genome mm10 \
--input_type h5 \
--output_type tsv
