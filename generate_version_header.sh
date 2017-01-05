#!/bin/bash
TEMPORARY_FILE="/tmp/$RANDOM"
VERSION_TEMPLATE="\xef\xbb\xbf
#pragma once

namespace principia {
namespace base {

char const BuildDate[] = \"%%DATE%%\";
char const Version[] =
    u8\"%%VERSION%%\";

}  // namespace base
}  // namespace principia"

echo -e "$VERSION_TEMPLATE" | 
sed "s/%%DATE%%/`date -d $(git log -1 --format=%cd --date=iso-strict) -u +%Y-%m-%dT%H:%M:%SZ`/" |
sed "s/%%VERSION%%/`git describe --tags --always --dirty --abbrev=40 --long`/" > $TEMPORARY_FILE
if cmp -s base/version.generated.h $TEMPORARY_FILE
then
  echo "No change to git describe, leaving base/version.generated.h untouched"
else
  echo "Updating base/version.generated.h"
  mv $TEMPORARY_FILE base/version.generated.h
fi
