#!/bin/bash
cd ./build
rm -r CMakeCache.txt CMakeFiles
cmake . -DCMAKE_CXX_COMPILER=g++
make clean
make -j4
cd ..
