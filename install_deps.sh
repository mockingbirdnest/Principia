#!/bin/bash

echo "apt-get installing: clang git unzip wget libc++-dev binutils make automake libtool curl cmake"
sudo apt-get install clang git unzip wget libc++-dev binutils make automake libtool curl cmake

mkdir -p deps
pushd deps

git clone "https://github.com/google/protobuf.git" --depth 1 -b "v3.0.0-alpha-1"
pushd protobuf
git am "../../documentation/Setup Files/protobuf.patch"
./autogen.sh
./configure CC=clang CXX=clang++ CXXFLAGS='-fPIC -m64 -std=c++11 -stdlib=libc++ -O3 -g' LDFLAGS='-stdlib=libc++'
make -j 8

popd
git clone https://github.com/Norgg/glog
pushd glog
# patch -p 1 -i "../../documentation/Setup Files/glog.patch"; true
./configure CC=clang CXX=clang++ CXXFLAGS='-fPIC -m64 -std=c++11 -stdlib=libc++ -O3 -g' LDFLAGS='-stdlib=libc++'
make -j 8

popd
wget https://googlemock.googlecode.com/files/gmock-1.7.0.zip
unzip gmock-1.7.0.zip
pushd gmock-1.7.0
patch -p 1 -i "../../documentation/Setup Files/gmock.patch"; true
./configure CC=clang CXX=clang++ CXXFLAGS='-fPIC -m64 -std=c++11 -I../glog/src -stdlib=libc++ -O3 -g' LDFLAGS='-stdlib=libc++'
make -j 8

popd
git clone https://github.com/pleroy/benchmark
pushd benchmark
cmake .
make

popd
git clone "https://chromium.googlesource.com/chromium/src.git" chromium -n --depth 1 -b "40.0.2193.1"
# $GitPromptSettings.RepositoriesInWhichToDisableFileStatus += join-path  (gi -path .).FullName chromium
pushd chromium
git config core.sparsecheckout true
cp "../../documentation/Setup Files/chromium_sparse_checkout.txt" .git/info/sparse-checkout
git checkout
git am "../../documentation/Setup Files/chromium.patch"

