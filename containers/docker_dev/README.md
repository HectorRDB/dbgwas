# Developing DBGWAS with docker

1. Install docker;
2. Download, create and start the container as:

	`docker run --name dbgwas_dev -it -v <dbgwas_root_dir>:/dbgwas leandroishilima/dbgwas:dev_0.5.4`

	Where `<dbgwas_root_dir>` is the root of DBGWAS in your machine. This will put you already in the mounted dbgwas folder in the container, ready for you to work.

3. If you are on **Linux**:
	* You can proceed by doing the common compilation commands on the mounted dbgwas folder:

	`mkdir build && cd build`

	`cmake .. && make && cd DBGWAS && make package`

3. If you are on **Mac/Windows**:
	* You **SHOULD NOT** compile on the `/dbgwas` mounted folder as you would do in Linux, since this is very slow in Mac/Windows;
	* Build the project somewhere else, outside of the mounted folder (here we build it on `/dbgwasDockerBuild`):

	`mkdir /dbgwasDockerBuild && cd /dbgwasDockerBuild`

	`cmake /dbgwas && make && cd DBGWAS && make package`

	* If you want to see the results using your browser, then you can copy to the mounted folder (e.g. to `/dbgwas/output`)

4. Note that the first compilation will take a lot of time. The next ones will be only incremental;
5. Do not delete anything after finishing your work, just exit the container;
6. This container can also run DBGWAS (it has R and the required packages installed);
7. After you exited the container, to go back to it, just run:

	`docker start -ai dbgwas_dev`