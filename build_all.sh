#!/bin/bash

# Build dependencies of C++ code
./build_deps.sh

# Build C++ code
make

# Build C# code
./build_ksp_plugin_adapter.sh