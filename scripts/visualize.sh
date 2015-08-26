#!/bin/bash
for i in {0..4}
do
    for j in {0..4}
    do
    python ./scripts/animation.py ./results/pl1_paths_${i}_${j}.txt ./results/pl2_paths_${i}_${j}.txt
    mv simulation.mp4 simulation_${i}_${j}.mp4
    done
done

