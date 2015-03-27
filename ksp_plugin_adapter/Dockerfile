FROM ubuntu:14.10

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install monodevelop -y
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb http://download.mono-project.com/repo/debian wheezy main" | tee /etc/apt/sources.list.d/mono-xamarin.list
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install mono-mcs -y

ADD . /opt/principia_adapter
WORKDIR /opt/principia_adapter
CMD mdtool build ksp_plugin_adapter_linux.csproj
