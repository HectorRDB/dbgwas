# Invokes docker to build portable dbgwas.
# Inspired by Páll Melsted blog (https://pmelsted.wordpress.com/2015/10/14/building-binaries-for-bioinformatics/),
# on how he and the other authors managed to make kallisto (Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter,
# Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519)
# portable in different linux distributions.
set -eu
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DBGWAS_DIR="$(dirname "${SCRIPT_DIR}")"

cd $DBGWAS_DIR

if [ -d "${DBGWAS_DIR}/build" ]; then
  echo "Please remove ${DBGWAS_DIR}/build before proceeding."
  exit 1
fi

sudo docker run -t -i --rm \
  -v ${DBGWAS_DIR}:/io \
  phusion/holy-build-box-64:2.0.1 \
  bash /io/portable_binary_builder/build_portable_dbgwas_core.sh
