#!/bin/bash
export OMP_NUM_THREADS=8
./build/opt lists/TerranTest.txt lists/ProtossTest.txt 100 20 10 10 
sh ./scripts/visualize.sh
