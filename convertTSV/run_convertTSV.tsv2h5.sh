#!/usr/bin/env bash

python convertTSV.py \
--input ../data/outputs/matrix.tsv \
--output ../data/outputs/matrix.h5 \
--genome mm10 \
--input_type tsv \
--output_type h5
