
#!/bin/bash
cd ~/Repos/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
cmake . -DCMAKE_CXX_COMPILER=g++
make clean
make -j4
cd ..
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="dynamic,1"
./build/opt lists/ZergTest.txt lists/ProtossTest.txt 25 5 25 5 
#sh ./scripts/visualize.sh

