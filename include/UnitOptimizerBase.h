#ifndef _UNIT_OPTIMIZER_BASE_H
#define _UNIT_OPTIMIZER_BASE_H

#include<limits>
#include<string>
#include<list>
#include<algorithm>
#include<vector>
#include<cmath>
#include<future>
#include<vector>

#include "Unit.h"
#include "UnitGenes.h"
#include "Race.h"
#include "InitPlayerUnits.h"
#include "MicroSimulation.h"
#include "Utilities.h"

class UnitOptimizerBase
{
public:
    virtual void setSelectionRate(float selectionRate)=0;
    virtual void setReproductionRate(float reproductionRate)=0;
    virtual void setMutationRate(float mutationRate)=0;
    virtual float getSelectionRate() const=0;
    virtual float getReproductionRate() const=0;
    virtual float getMutationRate() const=0;
    virtual void resizePopulation1(size_t size)=0;
    virtual void resizePopulation2(size_t size)=0;
    virtual void clearPopulation1()=0;
    virtual void clearPopulation2()=0;
    virtual void printBest()=0;
    virtual void setFieldSize(double, double)=0;
    virtual void optimize(size_t iterations, size_t stepsPerIteration)=0;
    virtual void optimize(size_t iterations, size_t stepsPerIteration, size_t initalPopulationSize, float selectionRate, float reproductionRate, float mutationRate)=0;
    virtual pair<UnitGenes,UnitGenes> getOptimum()=0;
    virtual pair<UnitGenes,UnitGenes> getOldOptimum()=0;
    virtual pair<double,double> getHighscores() const = 0;
    virtual pair<vector<double>,vector<double>> getSurvivalAverage() const = 0;

};

#endif
