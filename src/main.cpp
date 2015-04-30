#include <iostream>
#include <string>
#include <future>
#include <QApplication>
#include <vector>
#include <random>
#include <chrono>


#include "../include/TemplateInit.h"

#include "../include/InitPlayerUnits.h"
#include "../include/UnitFactory.h"
#include "../include/PlayerState.h"
#include "../include/MicroSimulation.h"
#include "../include/UnitOptimizer.h"
#include "../include/UnitOptimizerBase.h"
#include "../include/raceSelector.h"
#include "../include/BaseSelector.h"
#include "../include/GuiInterface.h"

using std::vector;
using std::string;
using std::initializer_list;

typedef struct OutputParameters
{
    UnitGenes initialGenes;
    pair<double,double> scores;
    float selectionRate;
    float reproductionRate;
    float mutationRate;
    vector<double> survivalAverage1;
    vector<double> survivalAverage2;
    OutputParameters(UnitGenes const& init, pair<double,double> sc, float sel, float rep, float mut, pair<vector<double>, vector<double>> survivalAverage)
	: initialGenes(init), scores(sc), selectionRate(sel), reproductionRate(rep), mutationRate(mut), survivalAverage1(survivalAverage.first), survivalAverage2(survivalAverage.second)
	{}
	    
} OutputParameters;

// Base path
string const baseDir = "";

static BaseSelector *microSim;

void writeOutput(SimulationResult const &, SimulationResult const &, OutputParameters const&);

bool run(size_t const b)
{
    return microSim->run(b);
}

void parseInputFile(std::string file, std::vector<int> &pl1, std::vector<int> &pl2)
{
    std::ifstream input(file.c_str(), std::ifstream::in);
    std::string tmp;
    input >> tmp;
    input >> tmp;
    for (int i = 0; i < 5; ++i)
    {
        input >> tmp;
        pl1.push_back(std::atoll(tmp.c_str()));
    }
    input >> tmp;
    for (int i = 0; i < 7; ++i)
    {
        input >> tmp;
        pl1.push_back(std::atoll(tmp.c_str()));
    }
    input >> tmp;
    input >> tmp;
    input >> tmp;
    input >> tmp;
    input >> tmp;
    input >> tmp;
    input >> tmp;
    input >> tmp;
    input >> tmp;
    input >> tmp;
    for (int i = 0; i < 5; ++i)
    {
        input >> tmp;
        pl2.push_back(std::atoll(tmp.c_str()));
    }
    input >> tmp;
    for (int i = 0; i < 7; ++i)
    {
        input >> tmp;
        pl2.push_back(std::atoll(tmp.c_str()));
    }
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: opt inputFile[Player1] inputFile[Player2]" << std::endl;
        return 1;
    }
    std::string file1 = argv[1];
    std::string file2 = argv[2];

    std::vector<std::string> units1;
    std::vector<std::string> units2;

    minimizeInputFile(units1, file1);
    minimizeInputFile(units2, file2);
    size_t subMe = file1.find_last_of('/');
    std::string tmp1 = file1.substr(subMe+1, file1.find_last_of(".txt")-4-subMe);
    subMe = file2.find_last_of('/');
    std::string tmp2 = file2.substr(subMe+1, file2.find_last_of(".txt")-4-subMe);

    std::string race1, race2;

    if (findRace(units1, race1) == false || findRace(units2, race2) == false)
    {
        std::cerr << "Couldn't determine Race. Abort." << std::endl;
        return 2;
    }

    std::string filePath1 = baseDir+"./data/"+race1+"/stats.txt";
    std::string filePath2 = baseDir+"./data/"+race2+"/stats.txt";
    std::string folder = "mkdir -p " + baseDir+"../results/";
    if (system(folder.c_str()))
    {}
    folder += race1 + "_" + race2;
    if (system(folder.c_str()))
    {}
    std::string outPath = baseDir+"./results/" + race1 + "_" + race2 + "/" + tmp1 + "_" + tmp2 + ".txt";
    std::cout << outPath << std::endl;
    size_t const iterations = 5, stepsPerIteration = 1000, initialPopulationSize = 100;

    int const x = (MAX-MIN)/2;

    UnitGenes initialGenes(x);
    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(0.0,1.0);
    float const selectionRate = 0.3; //dist(gen);
    float const reproductionRate = 0.75;//dist(gen);
    float const mutationRate = 0.75;//dist(gen);
    double fieldSize = 150;
    Vec2D minPos(0.0);
    Vec2D maxPos(fieldSize);

    UnitOptimizerBase *UOB;
    if (race1 == "Terran" && race2 == "Protoss")
    {
        UOB = new UnitOptimizer<Terran,Protoss>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Terran, Protoss>(minPos, maxPos, filePath1, filePath2);
    }
    else if (race1 == "Terran" && race2 == "Zerg")
    {
        UOB = new UnitOptimizer<Terran, Zerg>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Terran, Zerg>(minPos, maxPos, filePath1, filePath2);
    }
    else if (race1 == "Protoss" && race2 == "Zerg")
    {
        UOB = new UnitOptimizer<Protoss, Zerg>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Protoss, Zerg>(minPos, maxPos, filePath1, filePath2);
    }
    else if (race1 == "Protoss" && race2 == "Terran")
    {
        UOB = new UnitOptimizer<Protoss, Terran>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Protoss, Terran>(minPos, maxPos, filePath1, filePath2);
    }
    else if (race1 == "Zerg" && race2 == "Terran")
    {
        UOB = new UnitOptimizer<Zerg, Terran>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Zerg, Terran>(minPos, maxPos, filePath1, filePath2);
    }
    else if (race1 == "Zerg" && race2 == "Protoss")
    {
        UOB = new UnitOptimizer<Zerg, Protoss>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Zerg, Protoss>(minPos, maxPos, filePath1, filePath2);
    }
    else if (race1 == "Zerg" && race2 == "Zerg")
    {
        UOB = new UnitOptimizer<Zerg, Zerg>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Zerg, Zerg>(minPos, maxPos, filePath1, filePath2);
    }
    else if (race1 == "Protoss" && race2 == "Protoss")
    {
        UOB = new UnitOptimizer<Protoss, Protoss>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Protoss, Protoss>(minPos, maxPos, filePath1, filePath2);
    }
    else if (race1 == "Terran" && race2 == "Terran")
    {
        UOB = new UnitOptimizer<Terran, Terran>(units1, units2, filePath1, filePath2, initialPopulationSize, initialGenes);
        microSim = new GuiInterface<Terran, Terran>(minPos, maxPos, filePath1, filePath2);
    }
    else
    {
        return 4;
    }

    std::ofstream outF(outPath.c_str(), std::ofstream::out);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    UOB->optimize(iterations, stepsPerIteration, initialPopulationSize, selectionRate, reproductionRate, mutationRate);
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Elapsed time: " << elapsed.count() << " milliseconds" <<  std::endl;
    auto a = UOB->getOptimum();
    outF << "Player1:" << std::endl;
    outF << a.first << std::endl;
    outF << "Player2:" << std::endl;
    outF << a.second << std::endl;
    outF.close();

    microSim->initBothPlayers(units1, units2);
    microSim->setPlayer1Genes(UOB->getOptimum().first);
    microSim->setPlayer2Genes(UOB->getOptimum().second);
    microSim->setTracking(false);

    int retval = 0;
    QApplication app(argc, argv);
    app.setStyle("plastique");
    GuiWindow window;
    window.show();
    window.setStartFunc(run);
    microSim->setGuiWindow(&window);
    retval = app.exec();
    free(microSim);
    free(UOB);
    return retval;

}

void writeOutput(SimulationResult const &res1, SimulationResult const &res2, OutputParameters const &param)
{
    if (system("mkdir -p results"))
    {}
    std::ofstream sensitivityResult(baseDir+"./results/sensitivityResult.txt", std::ofstream::out);
    std::ofstream survivalResult(baseDir+"./results/survivalResult.txt", std::ofstream::out);
    std::ofstream pathPl1(baseDir+"./results/paths_pl1.txt", std::ofstream::out);
    std::ofstream pathPl2(baseDir+"./results/paths_pl2.txt", std::ofstream::out);
    std::ofstream resPl1(baseDir+"./results/results_pl1.txt", std::ofstream::out);
    std::ofstream resPl2(baseDir+"./results/results_pl2.txt", std::ofstream::out);

    sensitivityResult << "Initial Genes:" << std::endl;
    sensitivityResult << param.initialGenes << std::endl;
    sensitivityResult << "Highest score for Player 1: " << param.scores.first << std::endl;
    sensitivityResult << "Highest score for Player 2: " << param.scores.second << std::endl << std::endl;
    sensitivityResult << "Selection rate: " << param.selectionRate << std::endl;
    sensitivityResult << "Reproduction rate: " << param.reproductionRate << std::endl;
    sensitivityResult << "Mutation rate: " << param.mutationRate << std::endl;

    survivalResult << "Survival Average" << std::endl;
    for(size_t i = 0; i < std::min(param.survivalAverage1.size(), param.survivalAverage2.size()); ++i)
    {
	survivalResult << param.survivalAverage1[i] << "\t" << param.survivalAverage2[i] << std::endl;
    }
    for (size_t i = 0; i < res1.paths.back().size(); ++i)
    {
        for (size_t j = 0; j < res1.paths.size(); ++j)
        {
            pathPl1 << res1.paths.at(j).at(i).x << "\t" << res1.paths.at(j).at(i).y << "\t";
        }
        pathPl1 << std::endl;
    }
    for (size_t i = res1.paths.back().size(); i < res1.paths.at(0).size(); ++i)
    {
        for (size_t j = 0; j < res1.paths.size(); ++j)
        {
            if (res1.paths.at(j).size() <= i)
            {
                pathPl1 << res1.paths.at(j).back().x << "\t" << res1.paths.at(j).back().y << "\t";
            }
            else
            {
                pathPl1 << res1.paths.at(j).at(i).x << "\t" << res1.paths.at(j).at(i).y << "\t";
            }
        }
        pathPl1 << std::endl;
    }
    for (size_t i = 0; i < res2.paths.back().size(); ++i)
    {
        for (size_t j = 0; j < res2.paths.size(); ++j)
        {
            pathPl2 << res2.paths.at(j).at(i).x << "\t" << res2.paths.at(j).at(i).y << "\t";
        }
        pathPl2 << std::endl;
    }
    for (size_t i = res2.paths.back().size(); i < res2.paths.at(0).size(); ++i)
    {
        for (size_t j = 0; j < res2.paths.size(); ++j)
        {
            if (res2.paths.at(j).size() <= i)
            {
                pathPl2 << res2.paths.at(j).back().x << "\t" << res2.paths.at(j).back().y << "\t";
            }
            else
            {
                pathPl2 << res2.paths.at(j).at(i).x << "\t" << res2.paths.at(j).at(i).y << "\t";
            }
        }
        pathPl2 << std::endl;
    }
    resPl1 << "Minerals:\t" << res1.minerals << std::endl;
    resPl1 << "Vespene Gas:\t" << res1.gas << std::endl;
    resPl1 << "Total Damage:\t" << res1.damage << std::endl;
    resPl1 << "Damage in %:\t" << res1.damagePercent << std::endl;
    resPl1 << "Survivors:\t" << res1.survivors << std::endl;
    resPl2 << "Minerals:\t" << res2.minerals << std::endl;
    resPl2 << "Vespene Gas:\t" << res2.gas << std::endl;
    resPl2 << "Total Damage:\t" << res2.damage << std::endl;
    resPl2 << "Damage in %:\t" << res2.damagePercent << std::endl;
    resPl2 << "Survivors:\t" << res2.survivors << std::endl;

    resPl1.close();
    resPl2.close();
    pathPl1.close();
    pathPl2.close();
    sensitivityResult.close();
    survivalResult.close();
}
