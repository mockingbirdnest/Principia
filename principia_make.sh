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

set +e
make -j ${PARALLELISM} \
  --keep-going \
  bin/${PRINCIPIA_PLATFORM}/benchmark \
  bin/${PRINCIPIA_PLATFORM}/nanobenchmark \
  ${TARGET}
MAKE_RESULT=$?
set -e

git add *.png

if git diff --quiet HEAD *.png; then
  echo same;
else
  git config user.email "enrico.dandolo@mockingbirdnest.com"
  git config user.name "Enrico Dandolo"
  git commit -m "Update goldens for ${AGENT_OS} ${PRINCIPIA_PLATFORM}"
  BRANCH_NAME="goldens-$(date '+%Y%m%dT%H%M%S')-${AGENT_OS}-${PRINCIPIA_PLATFORM}"
  git checkout -b ${BRANCH_NAME}
  git push --set-upstream "https://${GH_TOKEN}@github.com/enrico-dandolo/Principia.git" ${BRANCH_NAME}
  gh pr create --fill --head enrico-dandolo:${BRANCH_NAME}
fi;

if [[ "${MAKE_RESULT}" != 0 ]]; then
  exit "${MAKE_RESULT}"
fi;

if [[ "${AGENT_OS?}" == "Darwin" ]]; then
  # See https://github.com/actions/virtual-environments/issues/2619#issuecomment-788397841
  # for why this is needed.
  sudo /usr/sbin/purge
fi
make release
