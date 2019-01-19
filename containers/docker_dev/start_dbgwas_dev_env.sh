#!/bin/bash
#just a small shortcut to docker run

#get the script and dbgwas dir
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
dbgwasDir="$(dirname "$(dirname "$scriptDir")")"

#run the container
#this will start the container and put you in the dbgwas folder, ready for you to work
docker run --name dbgwas_dev --rm -it -v ${dbgwasDir}:/dbgwas leandroishilima/dbgwas:dev_0.5.4

#you can proceed by doing the common compilation commands:
#mkdir build && cd build && cmake .. && make && cd DBGWAS/DBGWAS && make package

#note that the first compilation will take a lot of time. The next ones will be only incremental
#do not delete anything after finishing your work, just exit the container
#this is made so that the next time you work, the compilation is incremental

#this container can also run DBGWAS (it has R and the required packages installed)
