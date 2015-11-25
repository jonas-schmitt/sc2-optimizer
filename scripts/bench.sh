#!/bin/bash
cd ~/Repos/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
module load gcc/5.1.0
cmake . -DCMAKE_CXX_COMPILER=g++
make clean
make -j4
cd ..
rm -f bench.out
for i in `seq 1 8`;
do
    echo "Threads: $i"
    let "s = $i * 20"
    echo "Population size: $s"
    OMP_NUM_THREADS=$i ./build/opt lists/TerranTest.txt lists/ProtossTest.txt $s 20 10 10 >> bench.out
done
