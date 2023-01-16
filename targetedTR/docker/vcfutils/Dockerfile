FROM ubuntu:20.04

# Install dependencies
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get -qq install -y \
	autoconf \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    libncurses5 libncurses5-dev \
    pkg-config \
    wget \
    zlib1g-dev

# Make sources directory
RUN mkdir /sources

# Install VCFTools
WORKDIR /sources
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
RUN tar -xzvf vcftools-0.1.16.tar.gz
WORKDIR vcftools-0.1.16
RUN ./autogen.sh && ./configure && make && make install

# Install htslib (with tabix)
WORKDIR /sources
RUN wget -O htslib-1.8.tar.bz2 https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
RUN tar -xjvf htslib-1.8.tar.bz2
WORKDIR htslib-1.8/
RUN ./configure
RUN make
RUN make install

# Install bcftools
WORKDIR /sources
RUN wget -O bcftools-1.9.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
RUN tar -xjvf bcftools-1.9.tar.bz2
WORKDIR bcftools-1.9
RUN ./configure
RUN make
RUN make install