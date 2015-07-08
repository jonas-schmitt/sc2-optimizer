#!/bin/bash
#PBS -N sc2-optimizer
#PBS -l nodes=1:ppn=40
#PBS -l walltime=2:00:00
#PBS -M jonas.schmitt@fau.de
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err
module load cmake
module load gcc/4.9.2
module load openmpi/1.8.3-gcc
cd ~/sc2-optimizer
./configure
make
./exec.sh
