#!/bin/bash
docker run -t -v $PWD:/opt/principia/src principia make -j 4
