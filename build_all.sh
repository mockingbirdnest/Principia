#!/bin/bash

# Build dependencies of C++ code
./build_deps.sh

# Build C++ and C# code
make -j8
