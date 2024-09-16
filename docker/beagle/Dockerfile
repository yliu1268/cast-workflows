FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    git-all wget \
    curl \
    make \
    cmake \
    gcc \
    autoconf \
    automake \
    libtool \
    subversion \
    pkg-config \
    openjdk-21-jdk \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    libncurses5 libncurses5-dev \
    zlib1g-dev


#default-jre \


# Download and install gcloud
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin


RUN wget "https://faculty.washington.edu/browning/beagle/beagle.27May24.118.jar" -O beagle.jar

# Install htslib 1.20
WORKDIR /sources
RUN wget -O htslib-1.20.tar.bz2 https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
RUN tar -xjvf htslib-1.20.tar.bz2
WORKDIR htslib-1.20/
RUN ./configure
RUN make
RUN make install





