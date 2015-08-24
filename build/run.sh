#!/bin/bash
cd ~/Repos/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
cmake . -DCMAKE_CXX_COMPILER=g++
#make clean
make -j4
cd ..
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="static"
nice -n 19 nohup ./build/opt lists/TerranTest.txt lists/ProtossTest.txt 50 10 10 10 > opt.out 2> opt.err < /dev/null &
