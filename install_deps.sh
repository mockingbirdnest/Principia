#!/bin/bash
set -euo pipefail

echo "Required prerequisites for build: build-essential clang libc++-dev libc++abi-dev monodevelop subversion git"
echo "Required runtime dependencies: libc++1"

#sudo apt-get install clang git unzip wget libc++-dev binutils make automake libtool curl cmake subversion

mkdir -p deps
pushd deps

for repo in protobuf glog googletest gipfeli abseil-cpp compatibility benchmark zfp; do
  if [ ! -d "$repo" ]; then
    git clone "https://github.com/mockingbirdnest/$repo.git"
  fi
  pushd "$repo"
  git checkout master
  git pull

  AGENT_OS=$(uname -s)

  PRINCIPIA_C_FLAGS="-fPIC -O3 -g -DNDEBUG"
  PRINCIPIA_CXX_FLAGS="-std=c++17"
  PRINCIPIA_LD_FLAGS="-stdlib=libc++"
  PRINCIPIA_MACOS_CXX_FLAGS="-D_LIBCPP_STD_VER=16"
  PRINCIPIA_MACOS_VERSION_MIN="10.12"

  if [ -f ./principia_variable_overrides.sh ]; then
    . ./principia_variable_overrides.sh
  fi

  if [ "${AGENT_OS?}" == "Darwin" ]; then
    C_FLAGS="${PRINCIPIA_C_FLAGS?} -mmacosx-version-min=${PRINCIPIA_MACOS_VERSION_MIN?} -arch x86_64"
    CXX_FLAGS="${PRINCIPIA_CXX_FLAGS?} ${PRINCIPIA_MACOS_CXX_FLAGS?}"
  elif [ "${AGENT_OS?}" == "Linux" ]; then
    C_FLAGS="${PRINCIPIA_C_FLAGS?} -m64"
    CXX_FLAGS="${PRINCIPIA_CXX_FLAGS?}"
  else
    C_FLAGS="${PRINCIPIA_C_FLAGS?}"
    CXX_FLAGS="${PRINCIPIA_CXX_FLAGS?}"
  fi
  LD_FLAGS="${C_FLAGS?} ${PRINCIPIA_LD_FLAGS?}"
  CXX_FLAGS="${CXX_FLAGS?} ${LD_FLAGS?}"

  if [ -f ./principia_make.sh ]; then
    . ./principia_make.sh
  fi

  popd
done

if [ ! -d "Optional" ]; then
  mkdir Optional
fi
pushd Optional
curl "https://raw.githubusercontent.com/llvm-mirror/libcxx/52f9ca28a39aa02a2e78fa0eb5aa927ad046487f/include/optional" > principia_optional_impl
touch __undef_macros
popd

popd
