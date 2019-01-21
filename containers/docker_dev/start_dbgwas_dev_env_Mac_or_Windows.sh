#!/bin/bash
# THIS IS LIKE start_dbgwas_dev_env_Linux.sh, BUT WITH TWO BIG MODIFICATIONS
# PARAMETER --rm IS REMOVED FROM THE DOCKER COMMAND, SO THE CONTAINER IS PERSISTED AFTER YOU EXIT
# RESTART THE CONTAINER AS docker start -ai dbgwas_dev AFTER YOU EXIT, KILL, STOP, ETC...
# THIS SHOULD BE USED ONLY ON MAC/WINDOWS MACHINES, AS DISK WRITING BETWEEN THESE OSes AND A LINUX CONTAINER IS SLOW
# ALSO, YOU ***SHOULD NOT*** COMPILE ON DBGWAS FOLDER, SINCE THIS IS A MOUNTED FOLDER AND THUS IT WILL BE VERY SLOW
# BUILD SOMEWHERE ELSE, OUTSIDE OF THE MOUNTED FOLDER (e.g. BUILD ON /dbgwasDockerBuild)
# IF YOU WANT TO SEE THE RESULTS USING YOUR BROWSER, THEN YOU CAN COPY TO THE MOUNTED FOLDER

#just a small shortcut to docker run

#get the script and dbgwas dir
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
dbgwasDir="$(dirname "$(dirname "$scriptDir")")"

#run the container
#this will download the container, if you don't have it, start it, put you in the dbgwas folder, ready for you to work
docker run --name dbgwas_dev -it -v ${dbgwasDir}:/dbgwas leandroishilima/dbgwas:dev_0.5.4

#you can proceed by doing the common compilation commands (NOTE THAT THIS COMPILE NOT IN THE MOUNTED FOLDER, AS THIS IS SLOW ON WINDOWS/MAC):
#mkdir /dbgwasDockerBuild && cd /dbgwasDockerBuild && cmake /dbgwas && make && cd DBGWAS/DBGWAS && make package

#note that the first compilation will take a lot of time. The next ones will be only incremental
#do not delete anything after finishing your work, just exit the container
#this is made so that the next time you work, the compilation is incremental

#this container can also run DBGWAS (it has R and the required packages installed)
# IF YOU WANT TO SEE THE RESULTS USING YOUR BROWSER, THEN YOU CAN COPY TO THE MOUNTED FOLDER