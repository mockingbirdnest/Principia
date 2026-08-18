#! /bin/bash
# It seems that protoc really wants its dependencies to be in /usr/local/lib.
# In some setups, e.g., Azure pipelines, this does not work, so we need to help
# it find its dynamic libraries.

readonly NEW_EXTENSION='.new'

if [[ "${PRINCIPIA_PLATFORM?}" != "x64" &&
      "${PRINCIPIA_PLATFORM?}" != "x64_AVX_FMA" ]]; then
  echo "PRINCIPIA_PLATFORM must be x64 or x64_AVX_FMA."
  exit 1
fi

if [[ "${PRINCIPIA_PLATFORM?}" == "x64_AVX_FMA" ]]; then
  platform_golden_suffix="_fma"
fi

if [[ "${AGENT_OS?}" == "Darwin" ]]; then
  os_golden_suffix="_macos"
  parallelism=$(sysctl -n hw.ncpu)
  target="each_test"
elif [[ "${AGENT_OS?}" == "Linux" ]]; then
  os_golden_suffix="_linux"
  export LD_LIBRARY_PATH="./deps/protobuf/src/.libs:$LD_LIBRARY_PATH"
  parallelism=$(nproc --all)
  target="each_package_test"
fi

echo "Parallelism is ${parallelism}."

make clean

# Build the target, catching errors and grabbing the status.
set +e
make -j ${parallelism} \
  --keep-going \
  bin/${PRINCIPIA_PLATFORM}/benchmark \
  bin/${PRINCIPIA_PLATFORM}/nanobenchmark \
  ${target}
make_exit_status=$?
set -e

echo "Make finished with status ${make_exit_status}."

# Add all PNG files so new files are tracked.
git add *.png
env

if git diff --quiet HEAD; then
  echo "No files changed."
else
  # `git diff --name-status --no-renames` outputs something like:
  # A       path/to/new.file
  # D       path/to/deleted.file
  # M       path/to/modified.file
  git config core.quotePath false
  git diff --name-status --no-renames HEAD \
      | awk '/\.png/ {
               if ($1 == "D") {
                 print("(deleted) " $2)
               } else {
                 system("shasum -b " $2)
               }
             }' \
      > /tmp/golden_hashes.txt
  cat /tmp/golden_hashes.txt
  final_hash=$(shasum /tmp/golden_hashes.txt | awk '{print($1)}')
  branch_name="goldens-${final_hash}-${AGENT_OS}-${PRINCIPIA_PLATFORM}"
  echo "Looking for branch ${branch_name}."
  echo "git ls-remote --exit-code --heads https://github.com/enrico-dandolo/Principia.git refs/heads/${branch_name}"

  # `ls-remote` fails if the branch does not exist, so we need to handle errors.
  set +e
  git ls-remote --exit-code                                   \
      --heads https://github.com/enrico-dandolo/Principia.git \
      refs/heads/${branch_name}
  #echo "Done looking $?"
  ls_remote_exit_code=$?
  set -e

  echo "Exit code is ${ls_remote_exit_code}"
  if (( ${ls_remote_exit_code} == 2 )); then
    git config user.email "enrico.dandolo@mockingbirdnest.com"
    git config user.name "Enrico Dandolo"
    git checkout -b ${branch_name}
    golden_suffix="${os_golden_suffix}${platform_golden_suffix}"
    files_changed="$(git diff --name-only HEAD | grep '\.png')"
    echo "Files changed ${files_changed}"
    for file in ${files_changed}; do
      echo "File ${file} was changed"
      if [[ -f ${file} ]]; then
        echo "Moving to ${file}.new"
        mv "${file}" "${file}${NEW_EXTENSION}"
      fi
      cp $(echo "${file}" | sed "s/${golden_suffix}\.png/.png/") "${file}"
    done
    git reset
    git add *.png
    git commit \
        -m "Copy default goldens over ${AGENT_OS} ${PRINCIPIA_PLATFORM} goldens"
    echo "Files changed ${files_changed}"
    for file in ${files_changed}; do
      echo "File ${file} was changed, again"
      if [[ -f "${file}${NEW_EXTENSION}" ]]; then
        echo "Moving ${file}.new"
        mv "${file}${NEW_EXTENSION}" "${file}"
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
        ${branch_name}
    echo "Creating PR"
    gh pr create --fill --head enrico-dandolo:${branch_name} \
        --title "Update goldens for ${AGENT_OS} ${PRINCIPIA_PLATFORM}"
  else
    echo "Branch ${branch_name} already exists."
    echo "Merge the following PR to update these goldens:"
    gh pr list --head ${branch_name}
  fi
fi

if (( "${make_exit_status}" != 0 )); then
  exit "${make_exit_status}"
fi

if [[ "${AGENT_OS?}" == "Darwin" ]]; then
  # See https://github.com/actions/virtual-environments/issues/2619#issuecomment-788397841
  # for why this is needed.
  sudo /usr/sbin/purge
fi

make release
