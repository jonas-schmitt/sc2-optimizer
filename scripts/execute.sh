#!/bin/bash
export OMP_NUM_THREADS=4
./build/opt lists/Terran.txt lists/Protoss.txt 50 10 5 5 -stats ./results 
./scripts/visualize.sh
