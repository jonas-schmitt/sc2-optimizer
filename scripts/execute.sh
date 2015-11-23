#!/bin/bash
export OMP_NUM_THREADS=8
./build/opt lists/TerranTest.txt lists/ProtossTest.txt 100 30 10 10 -stats ./results 
./scripts/visualize.sh
