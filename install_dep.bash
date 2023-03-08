#!/bin/bash

git clone https://github.com/wckdouglas/diffexpr.git &&

conda config --add channels conda-forge &&

conda create -y -n diffexpre python=3.7 && 

conda activate diffexpre && 

conda install -y -c conda-forge tzlocal &&

conda install -y -c conda-forge rpy2 &&
conda install -y -c conda-forge biopython &&

conda install -y -c conda-forge reportlab &&

conda install -y -c conda-forge pytest-cov &&

conda install -y -c bioconda bioconductor-deseq2 &&

conda install -y -c conda-forge codecov &&

pip3 install --upgrade pip &&


conda install -y -c anaconda dash-bio &&

Rscript ./diffexpr/setup.R &&

python ./diffexpr/setup.py install &&

pip install contextvars goatools mygene packaging


