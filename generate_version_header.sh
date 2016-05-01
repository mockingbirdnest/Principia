#!/bin/bash
VERSION_TEMPLATE="\xef\xbb\xbf
#pragma once
namespace principia {
namespace base {

  char const kBuildDate[] = \"%%DATE%%\";
  char const kVersion[] =
      u8\"%%VERSION%%\";

}  // namespace base
}  // namespace principia"

echo -e "$VERSION_TEMPLATE" | 
sed "s/%%DATE%%/`date -u +%Y-%m-%dT%H:%M:%SZ`/" |
sed "s/%%VERSION%%/`git describe --tags --always --dirty --abbrev=40 --long`/" > base/version.hpp
