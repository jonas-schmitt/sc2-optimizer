#!/bin/bash
export OMP_NUM_THREADS=8
mkdir -p results8
./build/opt lists/Terran.txt lists/Protoss.txt 10000 50 10 10 -stats ./results8
#./scripts/visualize.sh
