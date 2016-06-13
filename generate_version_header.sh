#!/bin/bash
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
sed "s/%%DATE%%/`date -u +%Y-%m-%dT%H:%M:%SZ`/" |
sed "s/%%VERSION%%/`git describe --tags --always --dirty --abbrev=40 --long`/" > base/version.generated.h
