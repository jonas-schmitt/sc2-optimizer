
#!/bin/bash
cd ~/Repos/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
cmake . -DCMAKE_CXX_COMPILER=g++
make clean
make -j4
cd ..
export OMP_NUM_THREADS=4
export OMP_SCHEDULE="dynamic,1"
./build/opt lists/TerranTest.txt lists/ProtossTest.txt 20 3 50 3 
#sh ./scripts/visualize.sh

