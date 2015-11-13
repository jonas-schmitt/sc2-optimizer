#!/bin/bash
export OMP_NUM_THREADS=4
./build/opt lists/TerranTest.txt lists/ProtossTest.txt 20 5 2 1 
#./scripts/visualize.sh
