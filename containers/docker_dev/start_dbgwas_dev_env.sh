#!/bin/bash

#get the dbgwas dir
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
dbgwasDir="$(dirname "$(dirname "$scriptDir")")"

#pull the image - we fixed at version 2.0.1 so that we are always stable
docker pull phusion/holy-build-box-64:2.0.1

#run the container
#this will start the container and put you in the dbgwas folder, ready for you to work
docker run --rm -it -v ${dbgwasDir}:/dbgwas phusion/holy-build-box-64:2.0.1 /bin/bash -c "cd /dbgwas && source /hbb_exe/activate && /bin/bash"

#you can proceed by doing the common compilation commands:
#mkdir build && cd build && cmake .. && make

#note that the first compilation will take a lot of time. The next ones will be only incremental
#do not delete anything after finishing your work, just exit the container
#this is made so that the next time you work, the compilation is incremental