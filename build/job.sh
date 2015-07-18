#!/bin/bash
#PBS -N sc2-optimizer
#PBS -l nodes=64:ppn=40
#PBS -l walltime=23:00:00
#PBS -M jonas.schmitt@fau.de
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err
module load gcc/4.9.2
module load intelmpi
cd ~/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
./configure.sh
make clean
make -j10
cd ..
export OMP_NUM_THREADS=40
/apps/rrze/bin/mpirun -pernode ./build/opt lists/TerranTest.txt lists/ProtossTest.txt 1024 100 100 100 > res64.out 
