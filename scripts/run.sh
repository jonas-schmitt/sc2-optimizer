#!/bin/bash
cd ~/Repos/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
module load gcc/5.1.0
cmake . -DCMAKE_CXX_COMPILER=g++-5
make clean
make -j4
cd ..
./scripts/execute.sh
