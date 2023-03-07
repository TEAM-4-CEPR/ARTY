Authors: ILANGO Guy 

Maintainer : ILANGO Guy

Version : 0.0.3  

Draw : Cezard Adeline

Affiliation: Research Center for Respiratory Diseases ; Team 4  (France)

Credit : CEPR ILANGO G., CEPR CEZARD A.



-----------------------------------------------
# Analyse youR daTa Yourself (ARTY)
-----------------------------------------------

ARTY is made for all kind of analysis. It just need a formated table to process it. 

Features : - PCA (2d , 3d ) 
           - Differential expression analysis (DE)
           - Heatmap
           - VolcanoPlot
           - Gene ontology
           
This tools can only be used with a GUI

------------------------------------------------
# Prerequisite
------------------------------------------------
- python 3.7
- tkinter
- numpy
- diffexpr
- goatools
- mygene
- plotly
- sklearn
- pandas
- rpy2

-------------------------------------------------
# Installation
-------------------------------------------------
## Linux 
install all package with install_dep.bash

## Windows 
Because diffexpr is not available on windows, you will need to use WSL2.
You can install it using powershell and type following command : 
```
# go into some folder into which you want the file to be downloaded
cd <somefolder>

# download Ubuntu 20.04
Invoke-WebRequest -Uri https://aka.ms/wslubuntu2004 -OutFile Ubuntu.appx -UseBasicParsing

# install downloaded *.appx file
Add-AppxPackage .\Ubuntu.appx
```
Then you have to install Xserver : https://sourceforge.net/projects/vcxsrv/
In WSL2 type following command : 
```
export DISPLAY=$(ip route|awk '/^default/{print $3}'):0.0
```

Then you have to install a linux conda with WSL2 : https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Finally run install_dep.bash           

-------------------------------------------------
# Diagram
-------------------------------------------------
![](./diagram.png)

-------------------------------------------------
# How to use
-------------------------------------------------

You have to load your count matrix with sample as column and gene as rows. You have to precise all replicate (ex  : A_1,A_2,A_3,B_1,B_2,B_3)
Then you can just click on the DE button, this will create automatically your design matrix. And output a DE.csv containing fold change and pvalue. 
The DE button use diffexpr which is a python implementation of DesEQ2. 

Volcano , PCA (2d , 3d ) : You have to load count matrix and design matrix , then just click on the button.

Heatmap , Gene ontology : you have to load the count matrix then just click on the button

