#!/bin/bash
#PBS -N sc2-optimizer
#PBS -l nodes=5:ppn=32
#PBS -l walltime=24:00:00
#PBS -q normal
#PBS -M jonas.schmitt@fau.de
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err
. /etc/profile.d/modules.sh
module load cmake
module load gcc/5.3
#module load intel/2016.0.109
module load  openmpi/1.10.2-gnu
#module load openmpi/1.10.2-intel
cd ~/Repos/sc2-optimizer/build
rm -rf CMakeCache.txt CMakeFiles
./configure.sh
make clean
make -j4
cd ..
mkdir -p ./results
OMP_NUM_THREADS=8 mpirun --npersocket 1 \
    -mca orte_num_sockets 4 -mca orte_num_cores 8 \
        ./build/opt lists/Zerg.txt lists/Protoss.txt 1000 50 10 10 -stats ./results > debug.out
