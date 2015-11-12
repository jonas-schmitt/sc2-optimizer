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

#include "../include/TemplateInit.h"
#include "../include/OptimizerInterface.h"
#include "../include/Race.h"


struct OptimizationParameter
{
    Vec2D minPos;
    Vec2D maxPos;
    std::string filePath1;
    std::string filePath2;
    size_t popSize;
    std::vector<std::string> const* buildOrder1;
    std::vector<std::string> const* buildOrder2;
    size_t nGoals;
    size_t tournamentSize;
    size_t crossover;
    size_t mutation;
    size_t iterations;
    size_t generations;
    int rank;
    int procs;
    int migrants;
};

template<typename Race1, typename Race2>
void runOptimization(OptimizationParameter const& p)
{
    for(int i = 0; i < 5; ++i)
    {
        bool const saveStatistics = i == 4 ? true : false;
        OptimizerInterface<Race1,Race2> opt(p.minPos, p.maxPos, p.filePath1, p.filePath2, p.popSize, *p.buildOrder1, *p.buildOrder2, p.nGoals);
        opt.optimize(p.tournamentSize, p.crossover, p.mutation, p.iterations, p.generations, p.rank, p.procs, p.migrants, saveStatistics);
        std::ofstream resFile("./results/res" + std::to_string(i) + ".dat");
        opt.determineWinner(resFile, p.rank, p.procs, saveStatistics);
        resFile.close();
    }
//    auto elapsed_min = std::chrono::duration_cast<std::chrono::minutes>(end - start);
//    auto elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);
//    if(p.rank == 0) std::cout << "Elapsed time: " << elapsed_min.count() << " min " << elapsed_sec.count() - elapsed_min.count()*60 << " sec" << std::endl;
}

int main(int argc, char *argv[])
{

    OptimizationParameter p;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &p.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p.procs);

    if(argc < 7)
    {
        if(p.rank == 0) std::cout << "Usage: opt buildOrder1 buildOrder2 population iterations generations goals" << std::endl;
        MPI_Finalize();
        return -1;
    }

    Terran terran;
    Zerg zerg;
    Protoss protoss;


    std::ifstream file1(argv[1]), file2(argv[2]);
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
        if(p.rank == 0) std::cout << "Error: Couldn't determine first Race" << std::endl;
        MPI_Finalize();
        return -1;
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
        if(p.rank == 0) std::cout << "Error: Couldn't determine first Race" << std::endl;
        MPI_Finalize();
        return -1;
    }
    p.filePath1 = "./data/"+race1;
    p.filePath2 = "./data/"+race2;
    p.buildOrder1 = &buildOrder1;
    p.buildOrder2 = &buildOrder2;


    p.popSize = atoi(argv[3]);
    p.iterations = atoi(argv[4]);
    p.generations = atoi(argv[5]);
    p.nGoals = std::min(p.popSize, static_cast<size_t>(atoi(argv[6])));
    p.migrants = std::max(static_cast<size_t>(10), p.popSize / 10);
    size_t constexpr minFieldSize = 100;
    size_t const fieldSize = std::max(minFieldSize, 10 * std::max(buildOrder1.size(), buildOrder2.size()));
    p.minPos = Vec2D(0.0);
    p.maxPos = Vec2D(fieldSize,fieldSize);


    if(race1 == "Terran" && race2 == "Terran")
    {

        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Terran,Terran>(p);
                }
            }
        }


    }
    else if(race1 == "Terran" && race2 == "Zerg")
    {

        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Terran,Zerg>(p);
                }
            }
        }
    }
    else if(race1 == "Terran" && race2 == "Protoss")
    {


        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Terran,Protoss>(p);
                }
            }
        }
    }
    else if(race1 == "Zerg" && race2 == "Terran")
    {

        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Zerg,Terran>(p);
                }
            }
        }
    }
    else if(race1 == "Zerg" && race2 == "Zerg")
    {

        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Zerg,Zerg>(p);
                }
            }
        }
    }
    else if(race1 == "Zerg" && race2 == "Protoss")
    {

        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Zerg,Protoss>(p);
                }
            }
        }
    }
    else if(race1 == "Protoss" && race2 == "Terran")
    {


        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Protoss,Terran>(p);
                }
            }
        }
    }
    else if(race1 == "Protoss" && race2 == "Zerg")
    {

        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Protoss,Zerg>(p);
                }
            }
        }
    }
    else if(race1 == "Protoss" && race2 == "Protoss")
    {


        for(size_t i = 2; i < 3; i+=2)
        {
            for(size_t j = 0; j < 1; ++j)
            {
                for(size_t k = 0; k < 1; ++k)
                {
                    p.tournamentSize = i;
                    p.crossover = j;
                    p.mutation = k;
                    runOptimization<Protoss,Protoss>(p);
                }
            }
        }
    }
    MPI_Finalize();
}

