FROM ubuntu:14.10

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y clang git zip unzip wget libc++-dev binutils make automake libtool subversion cmake curl

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb http://download.mono-project.com/repo/debian wheezy main" | tee /etc/apt/sources.list.d/mono-xamarin.list
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install monodevelop -y

RUN git config --global user.email "docker@example.com"
RUN git config --global user.name docker

ADD documentation/ /opt/principia/documentation/

WORKDIR /opt/principia/
ADD install_deps.sh /opt/principia/install_deps.sh
RUN ./install_deps.sh

ENV ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5

WORKDIR /opt/principia/src
ADD ksp_plugin_adapter/ /opt/principia/src/ksp_plugin_adapter/
RUN mkdir "/opt/principia/KSP Assemblies"
RUN cp ksp_plugin_adapter/*.dll "/opt/principia/KSP Assemblies/"

CMD make -j4 DEP_DIR=../deps
