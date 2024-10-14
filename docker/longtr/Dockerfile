FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    autoconf \
    curl \
    zlib1g-dev libbz2-dev libz-dev liblzma-dev bzip2 \
    libcurl4-openssl-dev libssl-dev \
    python3 python3-pip \
    make gcc g++ \
    unzip \
    libfuse-dev fuse \
    git-all wget \
    cmake

RUN git clone https://github.com/gymrek-lab/longtr
WORKDIR longtr/lib
RUN git clone https://github.com/rvaser/spoa
WORKDIR ..

RUN make
RUN cp LongTR /usr/bin/