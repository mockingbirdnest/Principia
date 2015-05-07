FROM ubuntu:14.10

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y clang git unzip wget libc++-dev binutils make automake libtool subversion cmake curl

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb http://download.mono-project.com/repo/debian wheezy main" | tee /etc/apt/sources.list.d/mono-xamarin.list
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install monodevelop -y

RUN git config --global user.email "docker@example.com"
RUN git config --global user.name docker

ADD documentation/ /opt/principia/documentation/

WORKDIR /opt/principia/
RUN git clone "https://github.com/google/protobuf.git" --depth 1 -b "v3.0.0-alpha-1"
WORKDIR /opt/principia/protobuf
RUN git am "../documentation/Setup Files/protobuf.patch"
RUN ./autogen.sh
RUN ./configure CC=clang CXX=clang++ CXXFLAGS='-fPIC -m64 -std=c++11 -stdlib=libc++ -O3 -g' LDFLAGS='-stdlib=libc++'
RUN make -j 8

WORKDIR /opt/principia/
RUN git clone https://github.com/Norgg/glog
WORKDIR /opt/principia/glog
RUN ./configure CC=clang CXX=clang++ CXXFLAGS='-fPIC -m64 -std=c++11 -stdlib=libc++ -O3 -g' LDFLAGS='-stdlib=libc++'
RUN make -j 8

WORKDIR /opt/principia/
RUN svn checkout http://googlemock.googlecode.com/svn/trunk/ gmock
RUN svn checkout http://googletest.googlecode.com/svn/trunk/ gtest

WORKDIR /opt/principia/gtest/
RUN wget "https://gist.githubusercontent.com/Norgg/241ee11d278c0a55cc96/raw/4b23a866c6631ba0077229be366e67cde18fb035/gtest_linux_thread_count.patch" -O thread_count.patch
RUN patch -p 0 -i thread_count.patch

WORKDIR /opt/principia/gmock
RUN patch -p 1 -i "../documentation/Setup Files/gmock.patch"; true

WORKDIR /opt/principia
RUN git clone https://github.com/pleroy/benchmark
WORKDIR /opt/principia/benchmark
RUN cmake .
RUN make

WORKDIR /opt/principia
RUN git clone "https://chromium.googlesource.com/chromium/src.git" chromium -n --depth 1 -b "40.0.2193.1"
WORKDIR /opt/principia/chromium
RUN git config core.sparsecheckout true
RUN cp "../documentation/Setup Files/chromium_sparse_checkout.txt" .git/info/sparse-checkout
RUN git checkout
RUN git am "../documentation/Setup Files/chromium.patch"

ENV ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5

WORKDIR /opt/principia
RUN mkdir "/opt/KSP Assemblies"
RUN cp ksp_plugin_adapter/*.dll "/opt/KSP Assemblies"

WORKDIR /opt/principia/src
CMD make -j4 DEP_DIR=..
