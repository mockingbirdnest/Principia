#! /bin/bash
# It seems that protoc really wants its dependencies to be in /usr/local/lib.
# In some setups, e.g., Azure pipelines, this does not work, so we need to help
# it find its dynamic libraries.
if [[ "${PRINCIPIA_PLATFORM?}" != "x64" &&
      "${PRINCIPIA_PLATFORM?}" != "x64_AVX_FMA" ]]; then
  echo "PRINCIPIA_PLATFORM must be x64 or x64_AVX_FMA."
  exit 1
fi

if [[ "${AGENT_OS?}" == "Darwin" ]]; then
  PARALLELISM=$(sysctl -n hw.ncpu)
  TARGET="each_test"
elif [[ "${AGENT_OS?}" == "Linux" ]]; then
  export LD_LIBRARY_PATH="./deps/protobuf/src/.libs:$LD_LIBRARY_PATH"
  PARALLELISM=$(nproc --all)
  TARGET="each_package_test"
fi

echo "Parallelism is ${PARALLELISM}."

make clean

make -j ${PARALLELISM} \
  bin/${PRINCIPIA_PLATFORM}/benchmark \
  bin/${PRINCIPIA_PLATFORM}/nanobenchmark \
  ${TARGET}

if [[ "${AGENT_OS?}" == "Darwin" ]]; then
  # See https://github.com/actions/virtual-environments/issues/2619#issuecomment-788397841
  # for why this is needed.
  sudo /usr/sbin/purge
fi
make release
