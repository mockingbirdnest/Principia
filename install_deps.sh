#!/bin/bash

echo "Required prerequisites for build: build-essential clang libc++-dev libc++abi-dev monodevelop subversion git"
echo "Required runtime dependencies: libc++1"

#sudo apt-get install clang git unzip wget libc++-dev binutils make automake libtool curl cmake subversion

BASE_FLAGS="-fPIC -O3 -g"
# determine platform for bitness

PLATFORM=$(uname -s)
if [ "$PLATFORM" == "Darwin" ]; then
    C_FLAGS="$BASE_FLAGS -mmacosx-version-min=10.7 -arch i386"
elif [ "$PLATFORM" == "Linux" ]; then
	BITNESS=$(uname -m)
	if [ "$BITNESS" == "x86_64" ]; then
	  C_FLAGS="$BASE_FLAGS -m64"
        else
      	   C_FLAGS="$BASE_FLAGS -m32"
        fi
else
    C_FLAGS="$BASE_FLAGS"
fi

LD_FLAGS="$C_FLAGS -stdlib=libc++"
CXX_FLAGS="-std=c++1y $LD_FLAGS"

mkdir -p deps
cd deps

git clone "https://github.com/mockingbirdnest/protobuf.git"
pushd protobuf
./autogen.sh
if [ "$PLATFORM" == "Linux" ]; then
    ./autogen.sh # Really definitely needs to run twice on Ubuntu for some reason.
fi
./configure CC=clang CXX=clang++ CXXFLAGS="$CXX_FLAGS" LDFLAGS="$LD_FLAGS" LIBS="-lc++ -lc++abi"
make -j8
popd

git clone "https://github.com/Norgg/glog"
pushd glog
./configure CC=clang CXX=clang++ CFLAGS="$C_FLAGS" CXXFLAGS="$CXX_FLAGS" LDFLAGS="$LD_FLAGS" LIBS="-lc++ -lc++abi"
make -j8
popd

git clone "https://github.com/mockingbirdnest/gmock"
pushd gmock
# gmock does not need to be compiled
popd
