#include <iostream>
#include <string>
#include <future>
#include <vector>
#include <random>
#include <chrono>


#include "../include/TemplateInit.h"
#include "../include/MicroSimulation.h"
#include "../include/soga.h"



int main(int argc, char *argv[])
{
    std::string filePath1 = "./data/Terran";
    std::string filePath2 = "./data/Zerg";
    Vec2D minPos(0.0,0.0);
    Vec2D maxPos(200.0, 200.0);
    MicroSimulation<Terran,Zerg> sim(minPos, maxPos, filePath1, filePath2);
}
