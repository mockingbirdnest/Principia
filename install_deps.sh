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
  PRINCIPIA_CXX_FLAGS="-std=c++17"
  PRINCIPIA_LD_FLAGS="-stdlib=libc++"
  # NOTE(egg): We need Clang 8, and therefore we use Xcode 11.  The libc++ provided by
  # Xcode 11 has the <filesystem> header, however its symbols are marked unavailable
  # unless macosx-version-min is at least 10.15 (Catalina).  Since we do not want to
  # require that version yet (quousque tandem?), we cannot use the libc++ filesystem.
  # Instead we inject our own (https://github.com/mockingbirdnest/Compatibility), and we
  # prevent libc++ from defining and relying on its own, by setting _LIBCPP_STD_VER.
  # However, if we set that to 14, we would be missing void_t and the variable templates
  # from <type_traits>, which are C++17 additions on which we heavily rely.
  # Luckily, the <type_traits> definitions are gated on _LIBCPP_STD_VER > 14, while the
  # <filesystem> usage is gated on _LIBCPP_STD_VER >= 17.  We can therefore get the
  # former without the latter by setting _LIBCPP_STD_VER to 16 ∈ ]14, 17[.  
  PRINCIPIA_MACOS_CXX_FLAGS="-D_LIBCPP_STD_VER=16"
  PRINCIPIA_MACOS_VERSION_MIN="10.12"
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

if [ ! -d "Optional" ]; then
  mkdir Optional
fi
pushd Optional
curl "https://raw.githubusercontent.com/llvm-mirror/libcxx/52f9ca28a39aa02a2e78fa0eb5aa927ad046487f/include/optional" > principia_optional_impl
touch __undef_macros
popd

popd
