FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    libjpeg8-dev zlib1g-dev \
    git

RUN git clone https://github.com/AI-sandbox/gnomix
WORKDIR gnomix
RUN git checkout 2eca7b0 # to make this reproducible

# Install requirements 
ENV SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
RUN DEBIAN_FRONTEND="noninteractive" apt-get install -qqy cmake
RUN pip3 install -r requirements.txt
RUN yes | pip3 uninstall numpy
RUN pip3 install numpy==1.24.4