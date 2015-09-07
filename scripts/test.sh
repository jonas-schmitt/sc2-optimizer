#!/bin/bash
cd ~/Repos/sc2-optimizer/build
make -j4
cd ..
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="dynamic,1"
./build/opt lists/TerranTest.txt lists/ProtossTest.txt 20 3 50 3 
#sh ./scripts/visualize.sh

