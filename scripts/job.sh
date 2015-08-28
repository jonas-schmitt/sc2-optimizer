#!/bin/bash
#PBS -N sc2-optimizer
#PBS -l nodes=4:ppn=40
#PBS -l walltime=23:00:00
#PBS -M jonas.schmitt@fau.de
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err
module load cmake
module load gcc/4.9.2
module load intelmpi
cd ~/sc2-optimizer/build
rm -r CMakeCache.txt CMakeFiles
./configure.sh
make clean
make -j10
cd ..
export OMP_NUM_THREADS=40
mkdir -p ./results
/apps/rrze/bin/mpirun -pernode ./build/opt lists/TerranTest.txt lists/ProtossTest.txt 100 10 100 10 > res.out 
