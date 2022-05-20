FROM mcr.microsoft.com/mirror/docker/library/ubuntu:20.04
LABEL author="Venkat S. Malladi"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y --fix-missing && \
  apt-get upgrade -y && \
  apt-get dist-upgrade -y && \
  apt-get autoremove -y && \
  apt-get update -y --fix-missing && \
  apt-get upgrade -y
RUN locale
RUN apt-get install -y gcc g++ apt-utils wget gzip pigz pbzip2 zip software-properties-common make


ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
  bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
  rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/miniconda/bin:${PATH}

RUN  conda config --add channels bioconda && \
  conda config --add channels conda-forge && \
    conda config --remove channels defaults && \
  conda update -n base -y conda


RUN conda install pangolin
RUN conda install bedtools
RUN conda install bowtie2
RUN conda install bwa
RUN conda install fastp=0.23.1
RUN conda install lofreq
RUN conda install bcftools
RUN conda install samtools
RUN conda install sambamba
RUN conda install trimmomatic
RUN conda install multiqc=1.12
RUN ln -s /miniconda/lib/libcrypto.so /miniconda/lib/libcrypto.so.1.0.0

WORKDIR /tmp/
