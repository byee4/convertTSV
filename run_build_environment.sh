#!/usr/bin/env bash

ENVNAME=convertTSV

conda create -y -n $ENVNAME python=2.7;
source activate $ENVNAME;
conda install -c conda-forge matplotlib=2.0.2 numpy=1.13.3 scipy=0.19.1 pandas=0.20.3 pytables=3.4.2;
