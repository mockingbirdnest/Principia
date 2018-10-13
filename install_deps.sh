#!/bin/bash
set -euo pipefail

echo "Required prerequisites for build: build-essential clang libc++-dev libc++abi-dev monodevelop subversion git"
echo "Required runtime dependencies: libc++1"

#sudo apt-get install clang git unzip wget libc++-dev binutils make automake libtool curl cmake subversion

BASE_FLAGS="-fPIC -O3 -g -DNDEBUG"
# determine platform for bitness

PLATFORM=$(uname -s)
if [ "$PLATFORM" == "Darwin" ]; then
    C_FLAGS="$BASE_FLAGS -mmacosx-version-min=10.11 -arch x86_64"
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
CXX_FLAGS="-std=c++14 $LD_FLAGS"

mkdir -p deps
cd deps

if [ ! -d "protobuf" ]; then
  git clone "https://github.com/mockingbirdnest/protobuf"
fi
pushd protobuf
git checkout master
git pull
./autogen.sh
if [ "$PLATFORM" == "Linux" ]; then
    ./autogen.sh # Really definitely needs to run twice on Ubuntu for some reason.
fi
./configure CC=clang CXX=clang++ CXXFLAGS="$CXX_FLAGS" LDFLAGS="$LD_FLAGS" LIBS="-lc++ -lc++abi"
make -j8
popd

if [ ! -d "glog" ]; then
  git clone "https://github.com/mockingbirdnest/glog"
fi
pushd glog
git checkout master
git pull
./autogen.sh
./configure CC=clang CXX=clang++ CFLAGS="$C_FLAGS" CXXFLAGS="$CXX_FLAGS" LDFLAGS="$LD_FLAGS" LIBS="-lc++ -lc++abi"
make -j8
popd

# googlemock/googletest don't need to be compiled
if [ ! -d "googlemock" ]; then
  git clone "https://github.com/mockingbirdnest/googlemock"
fi
pushd googlemock
git checkout master
git pull
popd

if [ ! -d "googletest" ]; then
  git clone "https://github.com/mockingbirdnest/googletest"
fi
pushd googletest
git checkout master
git pull
popd

if [ ! -d "gipfeli" ]; then
  git clone "https://github.com/mockingbirdnest/gipfeli"
fi
pushd gipfeli
git checkout master
git pull
make
popd

if [ ! -d "compatibility" ]; then
  git clone "https://github.com/mockingbirdnest/compatibility"
fi
pushd compatibility
git checkout master
git pull
popd

if [ ! -d "Optional" ]; then
  mkdir Optional
fi
pushd Optional
curl "https://raw.githubusercontent.com/llvm-mirror/libcxx/52f9ca28a39aa02a2e78fa0eb5aa927ad046487f/include/optional" > principia_optional_impl
touch __undef_macros
popd

if [ ! -d "benchmark" ]; then
  git clone "https://github.com/mockingbirdnest/benchmark"
fi
pushd benchmark
git checkout master
git pull
cmake -DCMAKE_C_COMPILER:FILEPATH=`which clang` -DCMAKE_CXX_COMPILER:FILEPATH=`which clang++` -DCMAKE_C_FLAGS="${C_FLAGS}" -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DCMAKE_LD_FLAGS="${LD_FLAGS}"
make -j8
popd
