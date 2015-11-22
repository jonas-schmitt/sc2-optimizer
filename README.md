# StarCraft II Optimizer
Building requirements: gcc >= 5.1.0, OpenMP, MPI
Additional requirements for the visualization: python-matplotlib, ffmpeg

To build the project, the following script can be used:
./scripts/build.sh

The binary can then be executed in the following way:
./bin/opt buildOrder1 buildOrder2 p m n s [-stats] [directory]

buildOrder1: Path to the file containing the first build order
buildOrder2: Path to the file containing the second build order
p: The size of the population of both players
The optimization is run for m x n generations
s: Number of strategies used for the objective function evaluation
After each nth generation the s best individuals are exchanged between both populations
-stats: If this argument is supplied, additional information about the optimization is gathered
directory: Path to the directory where the additional files should be saved, can be only specified in case of -stats. Default value: "./" 

In the directory "lists" a number of sample build orders are provided.
To enable upgrades for specific units, the respective configuration files in the directory data/Race/upgrades must be adapted.
"0" denotes that a specific upgrade is disabled, "1" that it is enabled.

In case -stats was specified, the stored unit paths can be visualized using the following python script:
python ./scripts/animation.py dir/pl1_paths_x_y.dat dir/pl2_paths_x_y.dat
where directory is the directory where the paths have been saved, as explained above



