#ifndef UNITOPTIMIZER_H
#define UNITOPTIMIZER_H

#include<limits>
#include<string>
#include<list>
#include<algorithm>
#include<vector>
#include<cmath>
#include<future>
#include<vector>
#include<unordered_set>

#include "Unit.h"
#include "UnitGenes.h"
#include "Race.h"
#include "InitPlayerUnits.h"
#include "UnitOptimizerBase.h"
#include "MicroSimulation.h"
#include "Utilities.h"

using std::string;
using std::vector;
using std::future;
using std::vector;
using std::unordered_set;


template<class Race1, class Race2>
    class UnitOptimizer : public UnitOptimizerBase
{
private:

    vector<string> mUnitList1;
    vector<string> mUnitList2;
    string mFilePath1;
    string mFilePath2;
    Vec2D mMinPos = Vec2D(0,0);
    Vec2D mMaxPos = Vec2D(150,150);
    vector<UnitGenes> mPopulation1;
    vector<UnitGenes> mPopulation2;
    unordered_set<size_t> mPopControl1;
    unordered_set<size_t> mPopControl2;
    float mSelectionRate = 0.25;
    float mReproductionRate = 0.5;
    float mMutationRate = 0.5;
    size_t mInitialPopulationSize;
    size_t const NTHREADS = std::thread::hardware_concurrency();
    int const DELTA = 5;
    vector<MicroSimulation<Race1,Race2> > mSim1;
    vector<MicroSimulation<Race2,Race1> > mSim2;
    pair<UnitGenes,UnitGenes> mOldOptimum;
    vector<double> mSurvival1;
    vector<double> mSurvival2;
    bool mOptFlag = false;

    template<typename T, typename U>
        void evaluate(vector<UnitGenes>& population, vector<MicroSimulation<T,U>>& simVec)
    {
        vector<future<void> > futureVec(NTHREADS);
        vector<vector<Fitness>> threadRes(NTHREADS);
        for(size_t i = 0;i < NTHREADS; ++i)
        {
            threadRes[i] = std::move(vector<Fitness>(population.size()));
            futureVec[i] = std::async(std::launch::async,
                                      [] (MicroSimulation<T,U>* sim, vector<UnitGenes>* population, vector<Fitness>* threadRes)
                                      {
                                          for(size_t i = 0; i < population->size(); ++i)
                                          {
                                              sim->setPlayer1Genes((*population)[i]);
                                              (*threadRes)[i] = sim->run(true);
                                              /* Remove this line for a maximum workload */
                                              //std::this_thread::sleep_for(std::chrono::microseconds(750));
                                          }
                                      }, &simVec[i],&population,&threadRes[i]);
        }
        std::for_each(population.begin(), population.end(), [](UnitGenes& genes) {genes.fitness = Fitness();});
        for(size_t i = 0; i < NTHREADS; ++i)
        {
            futureVec[i].get();
            for(size_t j = 0; j < population.size(); ++j)
            {
                population[j].fitness = population[j].fitness + threadRes[i][j]/NTHREADS;

            }

        }
    }
    void initialize(size_t const N);
    void mutate(size_t const N);
    void crossover(size_t const N);
    void resetPopulation();


public:
UnitOptimizer(vector<string> unitList1, vector<string> unitList2, string filePath1, string filePath2, size_t initialPopulationSize)
    : mUnitList1(unitList1), mUnitList2(unitList2), mFilePath1(filePath1), mFilePath2(filePath2),
        mInitialPopulationSize(initialPopulationSize)
        {
            mPopulation1.reserve(initialPopulationSize);
            mPopulation2.reserve(initialPopulationSize);
            int const x = (MAX-MIN)/2;

            UnitGenes start(x);

            for(size_t i = 0; i < NTHREADS; ++i)
            {
                mSim1.emplace_back(mMinPos, mMaxPos, mFilePath1, mFilePath2);
                mSim2.emplace_back(mMinPos, mMaxPos, mFilePath2, mFilePath1);
            }
            for(size_t i = 0; i < NTHREADS; ++i)
            {
                mSim1[i].initBothPlayers(unitList1, unitList2);
                mSim1[i].setPlayer2Genes(start);
                mSim2[i].initBothPlayers(unitList2, unitList1);
                mSim2[i].setPlayer2Genes(start);
            }
        }
UnitOptimizer(vector<string> unitList1, vector<string> unitList2, string filePath1, string filePath2, size_t initialPopulationSize, float selectionRate, float reproductionRate, float mutationRate)
    : UnitOptimizer(unitList1, unitList2, filePath1, filePath2, initialPopulationSize)
    {
        if(selectionRate < EPS || selectionRate > 1-EPS
           || reproductionRate < EPS || reproductionRate > 1-EPS
           || mutationRate < EPS || mutationRate > 1-EPS)
        {
            throw std::invalid_argument("UnitOptimizer::UnitOptimizer(...): All rates must be values in the interval (0,1)");
        }

        mSelectionRate = selectionRate;
        mReproductionRate = reproductionRate;
        mMutationRate = mutationRate;
    }

UnitOptimizer(vector<string> unitList1, vector<string> unitList2, string filePath1, string filePath2, size_t initialPopulationSize, UnitGenes const & initialGenes)
    : mUnitList1(unitList1), mUnitList2(unitList2), mFilePath1(filePath1), mFilePath2(filePath2),
        mInitialPopulationSize(initialPopulationSize)
        {
            mPopulation1.reserve(initialPopulationSize);
            mPopulation2.reserve(initialPopulationSize);


            for(size_t i = 0; i < NTHREADS; ++i)
            {
                mSim1.emplace_back(mMinPos, mMaxPos, mFilePath1, mFilePath2);
                mSim2.emplace_back(mMinPos, mMaxPos, mFilePath2, mFilePath1);
            }
            for(size_t i = 0; i < NTHREADS; ++i)
            {
                mSim1[i].initBothPlayers(unitList1, unitList2);
                mSim1[i].setPlayer2Genes(initialGenes);
                mSim2[i].initBothPlayers(unitList2, unitList1);
                mSim2[i].setPlayer2Genes(initialGenes);
            }
        }

    void setSelectionRate(float selectionRate)
    {
        mSelectionRate = selectionRate;
    }
    void setReproductionRate(float reproductionRate)
    {
        mReproductionRate = reproductionRate;
    }
    void setMutationRate(float mutationRate)
    {
        mMutationRate = mutationRate;
    }
    float getSelectionRate() const
    {
        return mSelectionRate;
    }
    float getReproductionRate() const
    {
        return mReproductionRate;
    }
    float getMutationRate() const
    {
        return mMutationRate;
    }

    void resizePopulation1(size_t size)
    {
        mPopulation1.resize(size);
    }
    void resizePopulation2(size_t size)
    {
        mPopulation2.resize(size);
    }

    void clearPopulation1()
    {
        mPopulation1.clear();
    }

    void clearPopulation2()
    {
        mPopulation2.clear();
    }

    void setFieldSize(double x, double y)
    {
        mMaxPos.x = x;
        mMaxPos.y = y;
    }

    void optimize(size_t iterations, size_t stepsPerIteration);
    void optimize(size_t iterations, size_t stepsPerIteration, size_t initalPopulationSize, float selectionRate, float reproductionRate, float mutationRate);

    void printBest();
    void printNBest(size_t N);
    void printHighestCounts(size_t N);

    pair<UnitGenes,UnitGenes> getOptimum();
    pair<UnitGenes,UnitGenes> getOldOptimum();
    pair<vector<UnitGenes>,vector<UnitGenes>> getHighestCounts(size_t N);
    pair<double,double> getHighscores() const;
    pair<vector<double>,vector<double>> getSurvivalAverage() const;
};



#endif // UNITOPTIMIZER_H
