#!/bin/bash

# untar your R installation. Make sure you are using the right version!
#tar -xzf R402.tar.gz
# (optional) if you have a set of packages (created in Part 1), untar them also
#tar -xzf mSSL_simus_util_packages.tar.gz

mkdir ./res_$1-$2
# make sure the script will use your R installation, 
# and the working directory as its home location
#export PATH=$PWD/R/bin:$PATH
#export RHOME=$PWD/R
#export R_LIBS=$PWD/packages

# run your script
Rscript p20q30n100.R $1 $2 $3

tar -czvf p20q30n100res_$1-$2.tar.gz ./res_$1-$2

