FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    git-all wget \
    curl \
    default-jre \
    make \
    cmake \
    gcc \
    autoconf \
    automake \
    libtool \
    subversion \
    pkg-config



RUN wget -qO- "https://cmake.org/files/v3.17/cmake-3.17.0-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C /usr/local

# Download and install gcloud
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin


# Install BEAGLE
RUN apt-get update && apt-get install -y cmake build-essential autoconf automake libtool pkg-config
RUN git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git 

RUN cd beagle-lib \
    mkdir build \
    cd build \
    cmake -DBUILD_SSE=OFF -DBUILD=_CUDA=OFF .. \
    make install \
    ldconfig

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
CMD java -jar beagle
CMD beagle


