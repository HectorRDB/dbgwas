#!/usr/bin/env bash

apt-get update
apt-get install -y r-base wget libidn11
printf 'install.packages("ape", repos="http://cran.us.r-project.org")\ninstall.packages("phangorn", repos="http://cran.us.r-project.org")\ninstall.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos=NULL, type="source")' > R_deps
Rscript R_deps
rm R_deps
mkdir dbgwas
cd dbgwas
wget https://www.dropbox.com/s/s9oojqfl1kgi4l5/DBGWAS-0.5.4-Linux-precompiled.tar.gz?dl=1 -O DBGWAS-0.5.4-Linux-precompiled.tar.gz
tar -zxvf DBGWAS-0.5.4-Linux-precompiled.tar.gz