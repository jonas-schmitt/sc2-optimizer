#!/bin/bash
cd ~/Repos/sc2-optimizer
cmake . -DCMAKE_CXX_COMPILER=g++
make clean
make -j4
nice -n 19 nohup ./opt lists/TerranAttack.txt lists/ProtossAttack.txt > opt.out 2> opt.err < /dev/null &
