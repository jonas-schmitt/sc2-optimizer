#!/bin/bash
export OMP_NUM_THREADS=8
./build/opt lists/TerranTest.txt lists/ProtossTest.txt 50 20 10 10 
./scripts/visualize.sh
