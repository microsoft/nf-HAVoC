FROM continuumio/miniconda3:4.9.2-alpine
LABEL author="Venkat S. Malladi"

RUN apk update && apk add git bash

RUN  conda config --add channels bioconda && \
  conda config --add channels conda-forge && \
  conda update -n base -y conda

RUN conda install bedtools
RUN conda install bowtie2
RUN conda install bwa
RUN conda install fastp=0.23.1
RUN conda install lofreq
RUN conda install bcftools
RUN conda install samtools
RUN conda install sambamba
RUN conda install trimmomatic
RUN conda install pangolin=3.1.17
RUN ln -s /opt/conda/lib/libcrypto.so /opt/conda/lib/libcrypto.so.1.0.0

WORKDIR /tmp/
