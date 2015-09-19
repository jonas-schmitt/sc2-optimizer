#!/bin/bash
#PBS -N sc2-optimizer
#PBS -l nodes=4:ppn=32
#PBS -l walltime=6:00:00
#PBS -q normal
#PBS -M jonas.schmitt@fau.de
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err
. /etc/profile.d/modules.sh
module load cmake
module load gcc/5.1.0
module load openmpi/1.6.5-ib
cd ~/Repos/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
./configure.sh
make clean
make -j4
cd ..
mkdir -p ./results
#OMP_NUM_THREADS=8 mpirun --npersocket 1 \
#    -mca orte_num_sockets 4 -mca orte_num_cores 8 \
#        ./build/opt lists/TerranTest.txt lists/ProtossTest.txt 50 100 50 10 > res.out
OMP_NUM_THREADS=8 mpirun -np 16 -bysocket --bind-to-socket \
        ./build/opt lists/TerranTest.txt lists/ProtossTest.txt 100 10 50 10 > res.out
