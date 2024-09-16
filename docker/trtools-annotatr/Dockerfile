FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    git-all \
    python3 \
    python3-pip \
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

# Insatll TRTools version with annotatr
RUN pip3 install --upgrade pip
RUN git clone https://github.com/gymrek-lab/trtools.git 
WORKDIR trtools
RUN git checkout annotatr-less-strict-2
RUN pip install -e .

# Make sure we have tabix
# Install htslib 1.20
WORKDIR /sources
RUN wget -O htslib-1.20.tar.bz2 https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
RUN tar -xjvf htslib-1.20.tar.bz2
WORKDIR htslib-1.20/
RUN ./configure
RUN make
RUN make install

# Install bcftools 1.20
WORKDIR /sources
RUN wget -O bcftools-1.20.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
RUN tar -xjvf bcftools-1.20.tar.bz2
WORKDIR bcftools-1.20
RUN ./configure
RUN make
RUN make install


