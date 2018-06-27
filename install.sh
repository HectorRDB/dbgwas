SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH="${SOURCE_DIR}/seqan/util/cmake" -DSEQAN_INCLUDE_PATH="${SOURCE_DIR}/seqan/include" ..
make -j 4
make package
cd ..