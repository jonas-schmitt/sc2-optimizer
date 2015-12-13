#!/bin/bash
export OMP_NUM_THREADS=8
./build/opt lists/Terran.txt lists/Protoss.txt 2500 20 10 10 -stats ./results 
./scripts/visualize.sh
