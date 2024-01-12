FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    git-all \
    python3 \
    python3-pip

RUN pip3 install --upgrade pip
RUN git clone https://github.com/gymrek-lab/trtools.git 
WORKDIR trtools
RUN git checkout f07a99d
RUN pip install -e .