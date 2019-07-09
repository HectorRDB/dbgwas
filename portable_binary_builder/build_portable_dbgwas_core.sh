#!/bin/bash

# Builds portable dbgwas.
# Inspired by Páll Melsted blog (https://pmelsted.wordpress.com/2015/10/14/building-binaries-for-bioinformatics/),
# on how he and the other authors managed to make kallisto (Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter,
# Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519)
# portable in different linux distributions.

set -e

# Activate Holy Build Box environment.
source /hbb_exe/activate

set -eux

#install git
yum install git -y

# compile and install dbgwas
cd io
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd DBGWAS/DBGWAS
make package
cp DBGWAS*Linux-precompiled.tar.gz /io/
