#!/bin/bash
conda config --add channels conda-forge &&

conda create -y -n diffexpr python=3.7 && 
conda install -y -c conda-forge tzlocal &&

conda install -y -c conda-forge rpy2 &&
conda install -y -c conda-forge biopython &&

conda install -y -c conda-forge reportlab &&

conda install -y -c conda-forge pytest-cov &&

conda install -y -c bioconda bioconductor-deseq2 &&

conda install -y -c conda-forge codecov &&

pip3 install --upgrade pip &&

pip install goatools &&

pip install mygene &&

pip install plotly &&



pip install packaging &&

conda install  -y -c anaconda dash-bio &&


pip install contextvars
