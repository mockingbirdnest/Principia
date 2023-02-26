#!/bin/bash
TEMPORARY_FILE="/tmp/${RANDOM}"
VERSION_TEMPLATE="\xef\xbb\xbf
#include \"base/version.hpp\"

namespace principia {
namespace base {
namespace _version {
namespace internal {

char const BuildDate[] = \"%%DATE%%\";
char const Version[] =
    u8\"%%VERSION%%\";

}  // namespace internal
}  // namespace _version
}  // namespace base
}  // namespace principia"

PLATFORM=$(uname -s)
if [ "$PLATFORM" == "Darwin" ]; then
  echo -e "$VERSION_TEMPLATE" |
  sed "s/%%DATE%%/`date -j -u -f \"%Y-%m-%d %H:%M:%S %z\" -u +%Y-%m-%dT%H:%M:%SZ \"$(git log -1 --format=%cd --date=iso)\"`/" |
  sed "s/%%VERSION%%/`git describe --tags --always --dirty --abbrev=40 --long`/" > $TEMPORARY_FILE
else
  echo -e "$VERSION_TEMPLATE" |
  sed "s/%%DATE%%/`date -d $(git log -1 --format=%cd --date=iso-strict) -u +%Y-%m-%dT%H:%M:%SZ`/" |
  sed "s/%%VERSION%%/`git describe --tags --always --dirty --abbrev=40 --long`/" > $TEMPORARY_FILE
fi

if cmp -s base/version.generated.cc $TEMPORARY_FILE
then
  echo "No change to git describe, leaving base/version.generated.cc untouched"
else
  echo "Updating base/version.generated.cc"
  mv $TEMPORARY_FILE base/version.generated.cc
fi
