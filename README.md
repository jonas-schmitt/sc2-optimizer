# sc2-optimizer

## Building
Requirements: gcc >= 5.1.0 (other compilers are currently untested), OpenMP, MPI  
Additional requirements for the visualization: python-matplotlib, ffmpeg  

To build the project under linux, the following script can be used:  
./scripts/build.sh  
Building on other operating systems have not been tested so far.  

## Running
The binary can then be executed in the following way:  
./build/opt buildOrder1 buildOrder2 p m n s [-stats] [dir]  

* buildOrder1: Path to the file containing the first build order
* buildOrder2: Path to the file containing the second build order
* p: The size of the population of both players
* m, n: The optimization is run for m x n generations. After each nth generation the s best individuals are exchanged between both populations
* s: Number of strategies used for the objective function evaluation
* -stats: If this argument is supplied, additional information about the optimization is gathered
* dir: Path to the directory where the additional files should be saved, can be only specified in case of -stats. Default value: "./"

For execution with multiple processes a wrapper like mpirun must be used.  
Different example scripts for executing the binary on a single node (run.sh and execute.sh) and on a cluster using PBS (submit.sh and job.sh) can be found in the folder "scripts".

In the directory "lists" a number of sample build orders are provided.  
To enable upgrades for specific units, the respective configuration files in the directory "data/Race/upgrades" must be adapted. 0 denotes that a specific upgrade is disabled, 1 that it is enabled.  

## Visualization
In case -stats was specified, the stored unit paths can be visualized using the following python script:  
python ./scripts/animation.py dir/pl1_paths_x_y.dat dir/pl2_paths_x_y.dat  

where dir is the directory where the paths have been saved, as explained above  



