# To build: docker build . -t leandroishilima/dbgwas:0.5.4
# However, prefer to pull it from dockerhub
FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive
COPY ./install.sh ./
RUN chmod +x ./install.sh && ./install.sh
ENTRYPOINT ["/dbgwas/bin/DBGWAS"]
CMD ["-help"]
