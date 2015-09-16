#!/bin/bash
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="dynamic,1"
./build/opt lists/TerranTest.txt lists/ProtossTest.txt 100 100 50 10
sh ./scripts/visualize.sh
