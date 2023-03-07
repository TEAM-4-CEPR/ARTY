#!/bin/bash
conda config --add channels conda-forge &&

conda create -n diffexpr python=3.7 && 
conda install -c conda-forge tzlocal &&

conda install -c conda-forge rpy2 &&
conda install -c conda-forge biopython &&

conda install -c conda-forge reportlab &&

conda install -c conda-forge pytest-cov &&

conda install -c bioconda bioconductor-deseq2 &&

conda install -c conda-forge codecov &&

pip3 install --upgrade pip &&

pip install goatools &&

pip install mygene &&

pip install plotly &&



pip install packaging &&

conda install -c anaconda dash-bio &&


pip install contextvars