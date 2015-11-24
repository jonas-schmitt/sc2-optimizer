#include <mpi.h>
#include <iostream>
#include <string>
#include <future>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <thread>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "../include/TemplateInit.h"
#include "../include/OptimizerInterface.h"
#include "../include/Race.h"

void exit_with_error(std::string error_msg, int rank)
{
    if(rank == 0)
    {
        std::cerr << error_msg << "\n";
        std::cout << "Usage: opt buildOrder1 buildOrder2 population iterations generations goals [-stats] [directory]" << std::endl;
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
}

struct OptimizationParameter
{
    Vec2D minPos;
    Vec2D maxPos;
    std::string filePath1;
    std::string filePath2;
    std::string dirPath;
    size_t popSize;
    std::vector<std::string> const* buildOrder1;
    std::vector<std::string> const* buildOrder2;
    size_t nGoals; // number of strategies used for fitness evaluation
    size_t tournamentSize;
    size_t crossover;
    size_t mutation;
    size_t iterations;
    size_t generations;
    int rank;
    int procs;
    int migrants;
    bool saveStatistics = false;
};

template<typename Race1, typename Race2>
void runOptimization(OptimizationParameter const& p)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    // Initialization
    OptimizerInterface<Race1,Race2> opt(p.minPos, p.maxPos, p.filePath1, p.filePath2, p.popSize, *p.buildOrder1, *p.buildOrder2, p.nGoals, p.dirPath);

    // Optimization
    opt.optimize(p.tournamentSize, p.crossover, p.mutation, p.iterations, p.generations, p.rank, p.procs, p.migrants, p.saveStatistics);

    std::ofstream resFile(p.dirPath+"/res.dat");

    // Final comparison between both populations
    opt.determineWinner(resFile, p.rank, p.procs, p.saveStatistics);

    end = std::chrono::system_clock::now();
    auto elapsed_min = std::chrono::duration_cast<std::chrono::minutes>(end - start);
    auto elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    if(p.rank == 0) std::cout << "Elapsed time: " << elapsed_min.count() << " min " << elapsed_sec.count() - elapsed_min.count()*60 << " sec" << std::endl;
    resFile.close();
}


int main(int argc, char *argv[])
{
    OptimizationParameter p;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &p.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p.procs);


    // Parse the command line arguments
    if(argc < 7 || argc > 9)
    {
        exit_with_error("Error: Invalid number of command line arguments", p.rank);
    }

    if(argc >= 8)
    {
        std::string arg(argv[7]);
        if(arg != "-stats")
        {
            exit_with_error("Error: Invalid command line arguments", p.rank);
        }
        p.saveStatistics = true;
    }

    Terran terran;
    Zerg zerg;
    Protoss protoss;


    std::ifstream file1(argv[1]), file2(argv[2]);
    if(!file1 || !file2)
    {
        exit_with_error("Error: Could not open build order file", p.rank);
    }
    std::string buf;
    std::vector<std::string> buildOrder1, buildOrder2;
    while(std::getline(file1, buf))
    {
        buildOrder1.push_back(buf);
    }
    while(std::getline(file2, buf))
    {
        buildOrder2.push_back(buf);
    }


    // Ensure that both build orders only contain valid units
    std::string race1, race2;
    if(std::find(terran.nameList.begin(), terran.nameList.end(), buildOrder1[0]) != terran.nameList.end())
    {
        race1 = "Terran";
    }
    else if(std::find(zerg.nameList.begin(), zerg.nameList.end(), buildOrder1[0]) != zerg.nameList.end())
    {
        race1 = "Zerg";

    }
    else if(std::find(protoss.nameList.begin(), protoss.nameList.end(), buildOrder1[0]) != protoss.nameList.end())
    {
        race1 = "Protoss";
    }
    else
    {
        exit_with_error("Error: Could not determine first Race", p.rank);
    }

    for(size_t i = 1; i < buildOrder1.size(); ++i)
    {
        if(race1 == "Terran" && std::find(terran.nameList.begin(), terran.nameList.end(), buildOrder1[i]) == terran.nameList.end())
        {
            if(p.rank == 0) std::cerr << "The unit " << buildOrder1[i] << " is not a " << race1 << " unit!" << std::endl;
            exit_with_error("Error: Invalid build order for player 1", p.rank);
        }
        else if(race1 == "Zerg" && std::find(zerg.nameList.begin(), zerg.nameList.end(), buildOrder1[i]) == zerg.nameList.end())
        {
            if(p.rank == 0) std::cerr << "The unit " << buildOrder1[i] << " is not a " << race1 << " unit!" << std::endl;
            exit_with_error("Error: Invalid build order for player 1", p.rank);
        }
        else if(race1 == "Protoss" && std::find(protoss.nameList.begin(), protoss.nameList.end(), buildOrder1[i]) == protoss.nameList.end())
        {
            if(p.rank == 0) std::cerr << "The unit " << buildOrder1[i] << " is not a " << race1 << " unit!" << std::endl;
            exit_with_error("Error: Invalid build order for player 1", p.rank);
        }
    }

    if(std::find(terran.nameList.begin(), terran.nameList.end(), buildOrder2[0]) != terran.nameList.end())
    {
        race2 = "Terran";
    }
    else if(std::find(zerg.nameList.begin(), zerg.nameList.end(), buildOrder2[0]) != zerg.nameList.end())
    {
        race2 = "Zerg";

    }
    else if(std::find(protoss.nameList.begin(), terran.nameList.end(), buildOrder2[0]) != protoss.nameList.end())
    {
        race2 = "Protoss";
    }
    else
    {
        exit_with_error("Error: Could not determine second Race", p.rank);
    }

    for(size_t i = 1; i < buildOrder2.size(); ++i)
    {
        if(race2 == "Terran" && std::find(terran.nameList.begin(), terran.nameList.end(), buildOrder2[i]) == terran.nameList.end())
        {
            if(p.rank == 0) std::cerr << "The unit " << buildOrder2[i] << " is not a " << race2 << " unit!" << std::endl;
            exit_with_error("Error: Invalid build order for player 2", p.rank);
        }
        else if(race2 == "Zerg" && std::find(zerg.nameList.begin(), zerg.nameList.end(), buildOrder2[i]) == zerg.nameList.end())
        {
            if(p.rank == 0) std::cerr << "The unit " << buildOrder2[i] << " is not a " << race2 << " unit!" << std::endl;
            exit_with_error("Error: Invalid build order for player 2", p.rank);
        }
        else if(race2 == "Protoss" && std::find(protoss.nameList.begin(), protoss.nameList.end(), buildOrder2[i]) == protoss.nameList.end())
        {
            if(p.rank == 0) std::cerr << "The unit " << buildOrder2[i] << " is not a " << race2 << " unit!" << std::endl;
            exit_with_error("Error: Invalid build order for player 2", p.rank);
        }
    }

    // Paths to the directories containing the unit stats
    p.filePath1 = "./data/"+race1;
    p.filePath2 = "./data/"+race2;

    p.buildOrder1 = &buildOrder1;
    p.buildOrder2 = &buildOrder2;


    p.popSize = atoi(argv[3]);
    if(p.popSize <= 0)
    {
        exit_with_error("Error: The population size must be greater zero", p.rank);
    }
    p.iterations = atoi(argv[4]);
    p.generations = atoi(argv[5]);
    if(p.iterations <= 0 ||  p.generations <= 0)
    {
        exit_with_error("Error: The total number of generations must be greater zero", p.rank);
    }
    p.nGoals = std::min(p.popSize, static_cast<size_t>(atoi(argv[6])));
    if(p.nGoals <= 0)
    {
        exit_with_error("Error: The number of strategies used for fitness evaluation must be greater zero", p.rank);
    }
    p.migrants = std::max(static_cast<size_t>(10), p.popSize / p.procs);
    p.dirPath = "./";

    if(argc == 9)
    {
        struct stat info;
        if(stat(argv[8], &info ) != 0)
        {
            exit_with_error("Error: Invalid directory", p.rank);
        }
        else if(info.st_mode & S_IFDIR)
        {
            p.dirPath = std::string(argv[8]);
        }
        else
        {
            exit_with_error("Error: Invalid directory", p.rank);
        }
    }

    // Dynamically adapt the field size according to the number of units
    size_t constexpr minFieldSize = 100;
    size_t const fieldSize = std::max(minFieldSize, 10 * std::max(buildOrder1.size(), buildOrder2.size()));

    p.minPos = Vec2D(0.0);
    p.maxPos = Vec2D(fieldSize,fieldSize);

    // The default is binary tournament selection, self-adaptive SBX and gaussian mutation
    size_t i = 2, j = 0, k = 0;

    // Run the optimization
    if(race1 == "Terran" && race2 == "Terran")
    {
        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Terran,Terran>(p);
    }
    else if(race1 == "Terran" && race2 == "Zerg")
    {
        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Terran,Zerg>(p);

    }
    else if(race1 == "Terran" && race2 == "Protoss")
    {
        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Terran,Protoss>(p);
    }
    else if(race1 == "Zerg" && race2 == "Terran")
    {
        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Zerg,Terran>(p);
    }
    else if(race1 == "Zerg" && race2 == "Zerg")
    {
        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Zerg,Zerg>(p);
    }
    else if(race1 == "Zerg" && race2 == "Protoss")
    {
        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Zerg,Protoss>(p);
    }
    else if(race1 == "Protoss" && race2 == "Terran")
    {

        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Protoss,Terran>(p);
    }
    else if(race1 == "Protoss" && race2 == "Zerg")
    {
        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Protoss,Zerg>(p);
    }
    else if(race1 == "Protoss" && race2 == "Protoss")
    {

        p.tournamentSize = i;
        p.crossover = j;
        p.mutation = k;
        runOptimization<Protoss,Protoss>(p);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}

