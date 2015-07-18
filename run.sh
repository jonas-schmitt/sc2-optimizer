#!/bin/bash
cd ~/Repos/sc2-optimizer
rm -r CMakeCache.txt CMakeFiles
cmake . -DCMAKE_CXX_COMPILER=g++
make clean
make -j4
export OMP_NUM_THREADS=1
export OMP_SCHEDULE="static"
nice -n 19 nohup mpirun -np 1 ./opt lists/TerranTest.txt lists/ProtossTest.txt 40 10 10 5 > opt.out 2> opt.err < /dev/null &
