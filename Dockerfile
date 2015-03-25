FROM ubuntu:14.10

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y clang git unzip wget
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y binutils make 

WORKDIR /opt/principia/
RUN wget https://github.com/google/protobuf/releases/download/v2.6.1/protobuf-2.6.1.tar.bz2
RUN tar xf protobuf-2.6.1.tar.bz2
WORKDIR /opt/principia/protobuf-2.6.1
RUN ./configure
RUN make -j 8

WORKDIR /opt/principia/
RUN wget https://googlemock.googlecode.com/files/gmock-1.7.0.zip
RUN unzip gmock-1.7.0.zip
WORKDIR /opt/principia/gmock-1.7.0
RUN ./configure
RUN make -j 8

WORKDIR /opt/principia/
RUN git clone https://github.com/google/glog
WORKDIR /opt/principia/glog
RUN ./configure
RUN make -j 8

WORKDIR /opt/principia
RUN git clone https://github.com/pleroy/benchmark
# WORKDIR /opt/principia/benchmark
# RUN ./configure
# RUN make -j 8

ADD . /opt/principia
WORKDIR /opt/principia/
RUN make 
