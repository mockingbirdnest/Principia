#! /bin/bash
# It seems that protoc really wants its dependencies to be in /usr/local/lib.
# In some setups, e.g., Azure pipelines, this does not work, so we need to help
# it find its dynamic libraries.
if [[ "${PRINCIPIA_PLATFORM?}" != "x64" &&
      "${PRINCIPIA_PLATFORM?}" != "x64_AVX_FMA" ]]; then
  echo "PRINCIPIA_PLATFORM must be x64 or x64_AVX_FMA."
  exit 1
fi

if [[ "${PRINCIPIA_PLATFORM?}" != "x64_AVX_FMA" ]]; then
  PLATFORM_GOLDEN_SUFFIX="_fma"
fi

if [[ "${AGENT_OS?}" == "Darwin" ]]; then
  OS_GOLDEN_SUFFIX="_macos"
  PARALLELISM=$(sysctl -n hw.ncpu)
  TARGET="each_test"
elif [[ "${AGENT_OS?}" == "Linux" ]]; then
  OS_GOLDEN_SUFFIX="_linux"
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

if git diff --quiet HEAD *.png; then
  echo same;
else
  git diff --compact-summary HEAD |
      awk '/\.png/ { if ($2 == "(gone)") { print("(deleted) " $1) } else { system("sha1sum -b " $1) } }' \
      > golden_hashes.txt
  HASH=$(sha1sum golden_hashes.txt | awk '{print($1)}')
  BRANCH_NAME="goldens-${HASH}-${AGENT_OS}-${PRINCIPIA_PLATFORM}"
  git ls-remote --exit-code                                   \
      --heads https://github.com/enrico-dandolo/Principia.git \
      refs/heads/${BRANCH_NAME}
  if [[$? == 2]]; then
    git config user.email "enrico.dandolo@mockingbirdnest.com"
    git config user.name "Enrico Dandolo"
    git checkout -b ${BRANCH_NAME}
    GOLDEN_SUFFIX="${OS_GOLDEN_SUFFIX}${PLATFORM_GOLDEN_SUFFIX?}"
    FILES_CHANGED="$(git diff --compact-summary HEAD | awk '/\.png/ { $1 }')"
    for file in ${FILES_CHANGED}; do
      if [[ -f ${file} ]]; then
        mv "${file}" "${file}.new"
      fi
      cp "$(echo "${file}" | sed "s/${GOLDEN_SUFFIX}\.png/.png/")" "${file}"
    done
    git add *.png
    git commit \
        -m "Copy default goldens over ${AGENT_OS} ${PRINCIPIA_PLATFORM} goldens"
    for file in ${FILES_CHANGED}; do
      if [[ -f "${file}.new" ]]; then
        mv "${file}.new" "${file}"
      else
        git rm "${file}"
      fi
    done
    git add *.png
    git commit -m "Update goldens for ${AGENT_OS} ${PRINCIPIA_PLATFORM}"
    git push --set-upstream                                           \
        "https://${GH_TOKEN}@github.com/enrico-dandolo/Principia.git" \
        ${BRANCH_NAME}
    gh pr create --fill --head enrico-dandolo:${BRANCH_NAME}
  else
    echo "Branch ${BRANCH_NAME} already exists."
    echo "Merge the following PR to update these goldens:"
    gh pr list --head ${BRANCH_NAME}
  fi;
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
