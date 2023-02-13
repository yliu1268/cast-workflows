FROM gcr.io/ucsd-medicine-cast/vcfutils:latest

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip

RUN pip3 install --upgrade pip
RUN pip3 install trtools==5.0.1

COPY trviz.py /usr/bin/
RUN chmod +x /usr/bin/trviz.py