#!/bin/bash
docker run -t --rm -v $PWD:/opt/principia/src principia #make -j4 DEP_DIR=.. ksp_plugin_test/test
