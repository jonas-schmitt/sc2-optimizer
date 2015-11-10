#!/bin/bash
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="dynamic,1"
./build/opt lists/TerranTest.txt lists/ProtossTest.txt 20 10 10 5
sh ./scripts/visualize.sh
