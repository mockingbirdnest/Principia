#!/bin/bash
# Please specify the KSP_MANAGED_FOLDER variable if the location is different from the default in the script
# Typically this is something like: "${LOCATION_OF_GAMES}/Kerbal Space Program/KSP_Data/Managed/"
# This folder contains the KSP and unity C# assembly to build against

if [ -z ${KSP_MANAGED_FOLDER} ]; then
	KSP_MANAGED_FOLDER="${HOME}/.steam/steam/steamapps/common/Kerbal Space Program/KSP_Data/Managed/"
fi

# This test will fail if any of the following us true:
# The compiler doesn't recognize the C++ standard
# libc++fs cannot be found
# std::filesystem is not available
clang++ -x c++ -stdlib=libc++ -std=c++1z -lc++fs - &>/dev/null <<_TESTFILE
#include <filesystem>
#include <iostream>
int main() {
	std::cout << std::filesystem::path("/");
};
_TESTFILE
filesystem_support_test=$(expr "$?" == "0")

if [ "${filesystem_support_test}" -eq "1" ]
then
	AssemblySearchPaths="/usr/lib/mono/3.5-api/;/usr/lib/mono/2.0-api/;${KSP_MANAGED_FOLDER}" CXX_SUPPORTS_STD_FILESYSTEM=1 make adapter
else
	AssemblySearchPaths="/usr/lib/mono/3.5-api/;/usr/lib/mono/2.0-api/;${KSP_MANAGED_FOLDER}" make adapter
fi
