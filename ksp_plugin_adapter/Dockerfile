FROM ubuntu:14.10

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install monodevelop -y
ADD . /opt/principia_adapter
WORKDIR /opt/principia_adapter
# RUN mdtool build ksp_plugin_adapter_linux.csproj
