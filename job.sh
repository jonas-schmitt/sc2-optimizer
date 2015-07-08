#!/bin/bash
#PBS -N sc2-optimizer
#PBS -l nodes=8:ppn=40
#PBS -l walltime=12:00:00
#PBS -M jonas.schmitt@fau.de
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err
module load gcc/4.9.2
module load intelmpi
cd ~/sc2-optimizer
rm -r CMakeCache.txt CMakeFiles
./configure.sh
make clean
make -j10
export OMP_NUM_THREADS=20
mpirun_rrze  -np 16 -pinexpr S0:0-19_S1:0-19_S2:0-19_S3:0-19_S4:0-19_S5:0-19_S6:0-19_S7:0-19_S8:0-19_S9:0-19_S10:0-19_S11:0-19_S12:0-19_S13:0-19_S14:0-19_S15:0-19 ./opt lists/TerranAttack.txt lists/ProtossAttack.txt > opt.out
