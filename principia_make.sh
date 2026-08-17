#! /bin/bash
# It seems that protoc really wants its dependencies to be in /usr/local/lib.
# In some setups, e.g., Azure pipelines, this does not work, so we need to help
# it find its dynamic libraries.
if [[ "${PRINCIPIA_PLATFORM?}" != "x64" &&
      "${PRINCIPIA_PLATFORM?}" != "x64_AVX_FMA" ]]; then
  echo "PRINCIPIA_PLATFORM must be x64 or x64_AVX_FMA."
  exit 1
fi

if [[ "${PRINCIPIA_PLATFORM?}" == "x64_AVX_FMA" ]]; then
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
  astronomy/test
#  ${TARGET}
MAKE_RESULT=$?
set -e
echo "Make finished with status ${MAKE_RESULT}"

# Add all PNG files so new files are tracked.
git add *.png

if git diff --quiet HEAD; then
  echo "No files changed"
else
  # git diff --name-status --no-renames outputs something like
  # A       path/to/new.file
  # D       path/to/deleted.file
  # M       path/to/modified.file
  git config core.quotePath false
  git diff --name-status --no-renames HEAD |
      awk '/\.png/ { if ($1 == "D") { print("(deleted) " $2) } else { system("shasum -b " $2) } }' \
      > golden_hashes.txt
  cat golden_hashes.txt
  HASH=$(shasum golden_hashes.txt | awk '{print($1)}')
  BRANCH_NAME="goldens-${HASH}-${AGENT_OS}-${PRINCIPIA_PLATFORM}"
  echo "Looking for branch ${BRANCH_NAME}..."
  echo "git ls-remote --exit-code --heads https://github.com/enrico-dandolo/Principia.git refs/heads/${BRANCH_NAME}"
  set +e
  git ls-remote --exit-code                                   \
      --heads https://github.com/enrico-dandolo/Principia.git \
      refs/heads/${BRANCH_NAME}
  #echo "Done looking $?"
  code=$?
  set -e
  echo "Exit code is ${code}"
  if [[ ${code} == 2 ]]; then
    git config user.email "enrico.dandolo@mockingbirdnest.com"
    git config user.name "Enrico Dandolo"
    git checkout -b ${BRANCH_NAME}
    GOLDEN_SUFFIX="${OS_GOLDEN_SUFFIX}${PLATFORM_GOLDEN_SUFFIX}"
    FILES_CHANGED="$(git diff --name-only HEAD | grep '\.png')"
    echo "Files changed ${FILES_CHANGED}"
    for file in ${FILES_CHANGED}; do
      echo "File ${file} was changed"
      if [[ -f ${file} ]]; then
        echo "Moving to ${file}.new"
        mv "${file}" "${file}.new"
      fi
      cp $(echo "${file}" | sed "s/${GOLDEN_SUFFIX}\.png/.png/") "${file}"
    done
    git reset
    git add *.png
    git commit \
        -m "Copy default goldens over ${AGENT_OS} ${PRINCIPIA_PLATFORM} goldens"
    echo "Files changed ${FILES_CHANGED}"
    for file in ${FILES_CHANGED}; do
      echo "File ${file} was changed, again"
      if [[ -f "${file}.new" ]]; then
        echo "Moving ${file}.new"
        mv "${file}.new" "${file}"
      else
        echo "Removing ${file}"
        git rm "${file}"
      fi
    done
    echo "Adding"
    git add *.png
    echo "Commiting"
    git commit -m "Update goldens for ${AGENT_OS} ${PRINCIPIA_PLATFORM}"
    echo "Pushing"
    git push --set-upstream                                           \
        "https://${GH_TOKEN}@github.com/enrico-dandolo/Principia.git" \
        ${BRANCH_NAME}
    echo "Creating PR"
    gh pr create --fill --head enrico-dandolo:${BRANCH_NAME}
  else
    echo "Branch ${BRANCH_NAME} already exists."
    echo "Merge the following PR to update these goldens:"
    gh pr list --head ${BRANCH_NAME} \
        --title "Update goldens for ${AGENT_OS} ${PRINCIPIA_PLATFORM}"
  fi
fi

if [[ "${MAKE_RESULT}" != 0 ]]; then
  exit "${MAKE_RESULT}"
fi

if [[ "${AGENT_OS?}" == "Darwin" ]]; then
  # See https://github.com/actions/virtual-environments/issues/2619#issuecomment-788397841
  # for why this is needed.
  sudo /usr/sbin/purge
fi
make release
