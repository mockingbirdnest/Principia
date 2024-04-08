#!/bin/bash
set -euo pipefail

echo "Required prerequisites for build: build-essential clang libc++-dev libc++abi-dev monodevelop subversion git"
echo "Required runtime dependencies: libc++1"

#sudo apt-get install clang git unzip wget libc++-dev binutils make automake libtool curl cmake subversion

mkdir -p deps
pushd deps

for repo in abseil-cpp benchmark config gipfeli glog googletest multiprecision protobuf zfp; do
  if [ ! -d "$repo" ]; then
    git clone "https://github.com/mockingbirdnest/$repo.git"
  fi
  pushd "$repo"
  git checkout master
  git pull

  # Azure pipelines define this variable for us.
  AGENT_OS=$(uname -s)

  # Pipeline variables.
  # Any changes made to the variables defined in this section must be reflected
  # in the variable group
  # https://dev.azure.com/mockingbirdnest/Principia/_library?itemType=VariableGroups&view=VariableGroupView&variableGroupId=1&path=Principia
  # and vice-versa.
  # Note that the Principia Makefile doesn't get these variables; flags are duplicated
  # there.
  PRINCIPIA_C_FLAGS="-fPIC -O3 -g -DNDEBUG"
  PRINCIPIA_CXX_FLAGS="-std=c++20"
  PRINCIPIA_LD_FLAGS="-stdlib=libc++"
  PRINCIPIA_MACOS_CXX_FLAGS="-D_LIBCPP_STD_VER=20"
  PRINCIPIA_MACOS_VERSION_MIN="10.13"
  # End pipeline variables.

  # Task group Make.
  # Any changes to this section must be reflected in the script for the task group Make,
  # https://dev.azure.com/mockingbirdnest/Principia/_taskgroup/7ac7ecad-ff96-4796-9870-36aa93f5bacf,
  # and vice-versa.
  if [ -f ./principia_variable_overrides.sh ]; then
    . ./principia_variable_overrides.sh
  fi

  if [ "${AGENT_OS?}" == "Darwin" ]; then
    C_FLAGS="${PRINCIPIA_C_FLAGS?} -mmacosx-version-min=${PRINCIPIA_MACOS_VERSION_MIN?} -arch x86_64"
    CXX_FLAGS="${PRINCIPIA_CXX_FLAGS?} ${PRINCIPIA_MACOS_CXX_FLAGS?}"
    export OSX_DEPLOYMENT_TARGET="${PRINCIPIA_MACOS_VERSION_MIN?}"
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
  # End task group Make.

  popd
done

popd
