#!/bin/bash
VERSION_TEMPLATE="#pragma once
namespace principia {
namespace base {
char const kBuildDate[] = \"%%DATE%%\";
char const kVersion[] = 
    \"%%VERSION%%\";
}  // namespace base
}  // namespace principia"

echo -e "$VERSION_TEMPLATE" | 
sed "s/%%DATE%%/`date -u`/" |
sed "s/%%VERSION%%/`git describe --tags --always --dirty --abbrev=40 --long`/" > base/version.hpp
