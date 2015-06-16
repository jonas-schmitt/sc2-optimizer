#include <iostream>
#include <string>
#include <future>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>


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

    size_t popSize = 10;
    size_t iterations = 2;
    size_t genPerIt = 2;
    Vec2D minPos(0.0), maxPos(200.0,200.0);

    if(race1 == "Terran" && race2 == "Terran")
    {
        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<SOGA<Terran,Terran>,SOGA<Terran,Terran> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
                }
            }
        }

    }
    else if(race1 == "Terran" && race2 == "Zerg")
    {
        for(size_t i = 0; i < 3; ++i)
        {
            for(size_t j = 0; j < 5; ++j)
            {
                for(size_t k = 0; k < 3; ++k)
                {
                    Optimizer<SOGA<Terran,Zerg>,SOGA<Zerg,Terran> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
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
                    Optimizer<SOGA<Terran,Protoss>,SOGA<Protoss,Terran> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
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
                    Optimizer<SOGA<Zerg,Terran>,SOGA<Terran,Zerg> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
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
                    Optimizer<SOGA<Zerg,Zerg>,SOGA<Zerg,Zerg> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
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
                    Optimizer<SOGA<Zerg,Protoss>,SOGA<Protoss,Zerg> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
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
                    Optimizer<SOGA<Protoss,Terran>,SOGA<Terran,Protoss> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
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
                    Optimizer<SOGA<Protoss,Zerg>,SOGA<Zerg,Protoss> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
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
                    Optimizer<SOGA<Protoss,Protoss>,SOGA<Protoss,Protoss> > opt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2);
                    opt.optimize(i, j, k, iterations, genPerIt);
                }
            }
        }
    }


}
