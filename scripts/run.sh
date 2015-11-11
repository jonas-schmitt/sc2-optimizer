#!/bin/bash
cd ~/Repos/sc2-optimizer/build
#rm -r CMakeCache.txt CMakeFiles
module load gcc/5.1.0
#cmake . -DCMAKE_CXX_COMPILER=g++
#make clean
make -j4
cd ..
nohup sh ./scripts/execute.sh > test.out 2> test.err < /dev/null &
