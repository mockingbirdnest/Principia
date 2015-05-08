#!/bin/bash

echo "Required prerequisites for build: build-essential clang libc++-dev monodevelop subversion git"
echo "Required runtime dependencies: libc++1"

#sudo apt-get install clang git unzip wget libc++-dev binutils make automake libtool curl cmake subversion

BASE_FLAGS="-fPIC -O3 -g"
# determine platform for bitness

PLATFORM=$(uname -s)
if [ "$PLATFORM" == "Darwin" ]; then
    C_FLAGS="$BASE_FLAGS -mmacosx-version-min=10.7 -arch i386"
elif [ "$PLATFORM" == "Linux" ]; then
	BITNESS=$(uname -p)
	if [ "$BITNESS" == "x86_64"]; then
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

git clone "https://github.com/google/protobuf.git" --depth 1 -b "v3.0.0-alpha-1"
pushd protobuf
git am "../../documentation/Setup Files/protobuf.patch"
./autogen.sh
./configure CC=clang CXX=clang++ CXXFLAGS="$CXX_FLAGS" LDFLAGS="$LD_FLAGS"
make -j8
popd

git clone https://github.com/Norgg/glog
pushd glog
./configure CC=clang CXX=clang++ CFLAGS="$C_FLAGS" CXXFLAGS="$CXX_FLAGS" LDFLAGS="$LD_FLAGS"
make -j8
popd

svn checkout http://googletest.googlecode.com/svn/trunk/ gtest
pushd gtest
wget "https://gist.githubusercontent.com/Norgg/241ee11d278c0a55cc96/raw/4b23a866c6631ba0077229be366e67cde18fb035/gtest_linux_thread_count.patch" -O thread_count.patch
patch -p 0 -i thread_count.patch
# gtest does not need to be compiled
popd

svn checkout http://googlemock.googlecode.com/svn/trunk/ gmock
pushd gmock
patch -p 1 -i "../../documentation/Setup Files/gmock.patch"; true
# gmock does not need to be compiled
popd

# there's a bug in benchmark on 32-bit
#git clone https://github.com/google/benchmark
#pushd benchmark
#cmake -DCMAKE_CXX_FLAGS="$CXX_FLAGS" .
#make
#popd

#git clone "https://chromium.googlesource.com/chromium/src.git" chromium -n --depth 1 -b "40.0.2193.1"
#pushd chromium
#git config core.sparsecheckout true
#cp "../../documentation/Setup Files/chromium_sparse_checkout.txt" .git/info/sparse-checkout
#git checkout
#git am "../../documentation/Setup Files/chromium.patch"
