#!/bin/bash
#PBS -N ...sc2-optimizer
#PBS -l nodes=5:ppn=32
#PBS -l walltime=24:00:00
#PBS -q normal 
#PBS -M jonas.schmitt@fau.de
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err
. /etc/profile.d/modules.sh
module load cmake
module load gcc/5.1.0
module load openmpi/1.6.5-ib
cp -r ~/Repos/data ~/Repos/sc2-optimizer/
cd ~/Repos/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
./configure.sh
make clean
make -j4
cd ..
mkdir -p ./results
OMP_NUM_THREADS=8 mpirun --npersocket 1 \
    -mca orte_num_sockets 4 -mca orte_num_cores 8 \
        ./build/opt lists/TerranTest.txt lists/ProtossTest.txt 250 20 10 10 -stats ./results > debug.out
