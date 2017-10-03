#!/usr/bin/env bash

python convertTSV.py \
--input ../data/outputs/matrix.mtx \
--rows ../data/outputs/matrix.mtx.rows \
--columns ../data/outputs/matrix.mtx.columns \
--output ../data/outputs/matrix.csv \
--input_type mtx \
--output_type csv
