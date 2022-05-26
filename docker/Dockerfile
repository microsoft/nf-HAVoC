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

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh && \
  bash Miniconda3-py38_4.12.0-Linux-x86_64.sh -p /miniconda -b && \
  rm Miniconda3-py38_4.12.0-Linux-x86_64.sh

ENV PATH=/miniconda/bin:${PATH}

RUN  conda config --add channels bioconda && \
  conda config --add channels conda-forge && \
    conda config --remove channels defaults

RUN conda install mamba -n base -c conda-forge
RUN mamba install -c conda-forge biopython=1.74
RUN mamba install -c bioconda pangolin=4.0.2
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

# upgrade multiqc to development to support Pangolin 4.0.2 1.13_dev version
RUN apt-get update -y --fix-missing
RUN apt-get install -y git
RUN pip install --upgrade --force-reinstall git+https://github.com/ewels/MultiQC.git
RUN ln -s /miniconda/lib/libcrypto.so /miniconda/lib/libcrypto.so.1.0.0

WORKDIR /tmp/
