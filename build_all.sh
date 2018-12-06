#!/bin/bash

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
./build_deps.sh

# Build C++ and C# code
if [ "${filesystem_support_test}" -eq "1" ]
then
	CXX_SUPPORTS_STD_FILESYSTEM=1 make -j8
else
	make -j8
fi