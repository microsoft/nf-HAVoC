FROM continuumio/miniconda3:latest
LABEL author="Venkat S. Malladi"

RUN conda config --add channels defaults && \
  conda config --add channels bioconda && \
  conda config --add channels conda-forge && \
  conda update -n base -c defaults -y conda

RUN conda install bcftools
RUN conda install bedtools
RUN conda install bowtie2
RUN conda install bwa
RUN conda install fastp=0.23.1
RUN conda install samtools
RUN conda install lofreq
RUN conda install sambamba
RUN conda install trimmomatic
RUN conda install pangolin
RUN ln -s /opt/conda/lib/libcrypto.so /opt/conda/lib/libcrypto.so.1.0.0
