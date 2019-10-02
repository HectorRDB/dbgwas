# To build: docker build . -t leandroishilima/dbgwas:dev_0.5.4
# However, prefer to pull it from dockerhub
FROM phusion/holy-build-box-64:2.0.1
RUN yum -y install R && echo -e 'install.packages("ape", repos="http://cran.us.r-project.org")\ninstall.packages("phangorn", repos="http://cran.us.r-project.org")\ninstall.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos=NULL, type="source")' > R_deps && Rscript R_deps && rm R_deps
ENTRYPOINT ["/bin/bash", "-c", "source /hbb_exe/activate && cd /dbgwas && /bin/bash"]