# Invokes docker to build portable dbgwas.
# Inspired by Páll Melsted blog (https://pmelsted.wordpress.com/2015/10/14/building-binaries-for-bioinformatics/),
# on how he and the other authors managed to make kallisto (Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter,
# Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519)
# portable in different linux distributions.
set -eux
cd ..
rm -rf build
sudo docker run -t -i --rm \
  -v `pwd`:/io \
  phusion/holy-build-box-64:latest \
  bash /io/portable_binary_builder/build_portable_dbgwas_core.sh
