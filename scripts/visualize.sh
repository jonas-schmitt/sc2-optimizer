#!/bin/bash
for i in {0..9}
do
    for j in {0..9}
    do
    python ./scripts/animation.py ./results/pl1_paths_${i}_${j}.dat ./results/pl2_paths_${i}_${j}.dat
    mv simulation.mp4 ./results/simulation_${i}_${j}.mp4
    done
done

