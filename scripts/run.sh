#!/bin/bash
cd ~/Repos/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
cmake . -DCMAKE_CXX_COMPILER=g++
make clean
make -j4
cd ..
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="dynamic,1"
nohup nice -n 19 ./build/opt lists/TerranTest.txt lists/ProtossTest.txt 100 10 100 10 > res.out 2> test.err < /dev/null &
#sh ./scripts/visualize.sh

