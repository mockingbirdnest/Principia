#!/bin/bash
# Pass the option "--skip-deps" if you don't want to rebuild the dependencies.
#
# Please specify the KSP_MANAGED_FOLDER variable if the location is different from the default in the script
# Typically this is something like: "${LOCATION_OF_GAMES}/Kerbal Space Program/KSP_Data/Managed/"
# This folder contains the KSP and unity C# assembly to build against

build_deps=1

if [ -z ${KSP_MANAGED_FOLDER} ]; then
	KSP_MANAGED_FOLDER="${HOME}/.steam/steam/steamapps/common/Kerbal Space Program/KSP_Data/Managed/"
fi

# Read command line arguments
for i in "$@"
do
case $i in
	--skip-deps)
	build_deps=0
	shift
	;;
	*)
	echo "Unknown option ${i} passed"
	exit 1
	;;
esac
done

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

# Build dependencies of C++ code
if [ ${build_deps} -eq 1 ]
then
	./build_deps.sh
fi

# Build C++ and C# code
if [ "${filesystem_support_test}" -eq "1" ]
then
	AssemblySearchPaths="/usr/lib/mono/3.5-api/;/usr/lib/mono/2.0-api/;${KSP_MANAGED_FOLDER}" CXX_SUPPORTS_STD_FILESYSTEM=1 make -j8
else
	AssemblySearchPaths="/usr/lib/mono/3.5-api/;/usr/lib/mono/2.0-api/;${KSP_MANAGED_FOLDER}" make -j8
fi