FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    bcftools \
    curl \
    git-all wget \
    make 



# Download and install gcloud
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin


# Download HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2
RUN tar -xvf htslib-1.19.tar.bz2
WORKDIR htslib-1.19
RUN ./configure --enable-gcs --enable-s3 --enable-libcurl
RUN make -j && make install