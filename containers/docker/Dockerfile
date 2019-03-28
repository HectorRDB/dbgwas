# To build: docker build . -t leandroishilima/dbgwas:0.5.2
# However, prefer to pull it from dockerhub
FROM ubuntu:17.10
RUN apt-get update && apt-get install -y r-base wget && printf 'install.packages("ape", repos="http://cran.us.r-project.org")\ninstall.packages("phangorn", repos="http://cran.us.r-project.org")\ninstall.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos=NULL, type="source")' > R_deps && Rscript R_deps && rm R_deps && mkdir dbgwas && cd dbgwas && wget https://www.dropbox.com/s/gf6n4ibcakcyo5k/DBGWAS-0.5.2-Linux-precompiled.tar.gz?dl=1 -O DBGWAS-0.5.2-Linux-precompiled.tar.gz && tar -zxvf DBGWAS-0.5.2-Linux-precompiled.tar.gz
ENTRYPOINT ["/dbgwas/bin/DBGWAS"]
CMD ["-help"]