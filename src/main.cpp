#include <iostream>
#include <string>
#include <future>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <thread>
#include <mpi.h>
#include <omp.h>

#include "../include/TemplateInit.h"
#include "../include/MicroSimulation.h"
#include "../include/soga.h"
#include "../include/Optimizer.h"
#include "../include/Race.h"

using std::cout;
using std::string;
using std::ifstream;

int main(int argc, char *argv[])
{



    if(argc < 3)
    {
        std::cout << "Usage: opt buildOrder1 buildOrder2" << std::endl;
        return -1;
    }
    int rank, procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    omp_set_num_threads(omp_get_max_threads());
    //cout << "Number of Threads used: " << omp_get_max_threads() << std::endl;

    Terran terran;
    Zerg zerg;
    Protoss protoss;

    string filePath1, filePath2;

    ifstream file1(argv[1]), file2(argv[2]);
    string buf;
    vector<string> buildOrder1, buildOrder2;
    while(std::getline(file1, buf))
    {
        buildOrder1.push_back(buf);
    }
    while(std::getline(file2, buf))
    {
        buildOrder2.push_back(buf);
    }


    string race1, race2;

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
        cout << "Error: Couldn't determine first Race" << std::endl;
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
        cout << "Error: Couldn't determine second Race" << std::endl;
        return -1;
    }
    filePath1 = "./data/"+race1;
    filePath2 = "./data/"+race2;


    size_t popSize = 100;
    size_t iterations = 100;
    size_t genPerIt = 10;
    size_t nGoals = 10;
    size_t migrants = popSize / 2;
    Vec2D minPos(0.0), maxPos(200.0,200.0);

//    MicroSimulation<Terran, Protoss> sim(minPos, maxPos, filePath1, filePath2);
//    sim.initBothPlayers(buildOrder1, buildOrder2);
//    Chromosome initChrom1(sim.getPlayer1ChromosomeLength ());
//    Chromosome initChrom2(sim.getPlayer2ChromosomeLength ());
//    for(auto& gene : initChrom1)
//    {
//        gene.flip(gene.size()-1);
//    }
//    for(auto& gene : initChrom2)
//    {
//        gene.flip(gene.size()-1);
//    }
//    sim.setPlayer1Chromosome (initChrom1);
//    sim.setPlayer2Chromosome (initChrom2);
//    Fitness fitness = sim.run (false);
//    cout << "Score: " << fitness.score << endl;
//    cout << "Damage: " << fitness.damage*100 << endl;
//    cout << "Health: " << fitness.health*100 << endl;
//    cout << "Minerals alive: " << fitness.minerals_alive*100 << endl;
//    cout << "Gas alive: " << fitness.gas_alive*100 << endl;
//    cout << "Minerals killed: " << fitness.minerals_killed*100 << endl;
//    cout << "Gas killed: " << fitness.gas_killed*100 << endl;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    if(race1 == "Terran" && race2 == "Terran")
    {
        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Terran,Terran, Player::first>,MOGA<Terran,Terran, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                }
            }
        }
        //opt.determineWinner(std::cout);

    }
    else if(race1 == "Terran" && race2 == "Zerg")
    {
        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Terran,Zerg, Player::first>,MOGA<Zerg,Terran, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                }
            }
        }
    }
    else if(race1 == "Terran" && race2 == "Protoss")
    {

        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Terran,Protoss, Player::first>,MOGA<Protoss,Terran, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                }
            }
        }
    }
    else if(race1 == "Zerg" && race2 == "Terran")
    {
        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Zerg,Terran, Player::first>,MOGA<Terran,Zerg, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                }
            }
        }
    }
    else if(race1 == "Zerg" && race2 == "Zerg")
    {
        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Zerg,Zerg, Player::first>,MOGA<Zerg,Zerg, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                }
            }
        }
    }
    else if(race1 == "Zerg" && race2 == "Protoss")
    {
        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Zerg,Protoss, Player::first>,MOGA<Protoss,Zerg, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                }
            }
        }
    }
    else if(race1 == "Protoss" && race2 == "Terran")
    {

        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Protoss,Terran, Player::first>,MOGA<Terran,Protoss, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                }
            }
        }
    }
    else if(race1 == "Protoss" && race2 == "Zerg")
    {
        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Protoss,Zerg, Player::first>,MOGA<Zerg,Protoss, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                }
            }
        }
    }
    else if(race1 == "Protoss" && race2 == "Protoss")
    {

        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<MOGA<Protoss,Protoss, Player::first>,MOGA<Protoss,Protoss, Player::second> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals);
                    opt.optimize(k, j, i, iterations, genPerIt, rank, procs, migrants);
                    opt.determineWinner(std::cout, rank, procs);
                }
            }
        }
    }

    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Elapsed time: " << elapsed.count() << " milliseconds" <<  std::endl;
    MPI_Finalize();


}
