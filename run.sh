#!/bin/bash
cd ~/Repos/sc2-optimizer
rm -r CMakeCache.txt CMakeFiles
cmake . -DCMAKE_CXX_COMPILER=g++
make clean
make -j4
nice -n 19 nohup ./opt lists/ProtossTest.txt lists/ProtossTest.txt > opt.out 2> opt.err < /dev/null &
