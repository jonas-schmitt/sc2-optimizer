#!/bin/bash
export OMP_NUM_THREADS=8
mkdir -p results4
./build/opt lists/Terran.txt lists/Protoss.txt 1000 50 10 100 -stats ./results4
#./scripts/visualize.sh
