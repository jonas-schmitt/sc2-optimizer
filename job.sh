#!/bin/bash
#PBS -N sc2-optimizer
#PBS -l nodes=1:ppn=40
#PBS -l walltime=2:00:00
#PBS -M jonas.schmitt@fau.de
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err
module load gcc/4.9.2
module load openmpi
cd ~/sc2-optimizer
rm -r CMakeCache.txt CMakeFiles
./configure.sh
make clean
make -j4
sh exec.sh
