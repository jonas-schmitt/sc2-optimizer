#include <string>
#include<vector>
#include<algorithm>
#include<unordered_set>
#include<random>
#include<chrono>
#include<limits>
#include<future>
#include<iostream>
#include<functional>
#include<thread>
#include<cmath>
#include<cstdio>


#include "UnitOptimizer.h"

using std::string;
using std::vector;
using std::future;
using std::async;
using std::function;

template<class Race1, class Race2>
void UnitOptimizer<Race1, Race2>::initialize(size_t const N)
{
    for(size_t i = std::min(mPopulation1.size(),mPopulation2.size()); i < N; ++i)
    {
        UnitGenes individual;
        while(mPopControl1.count(individual.getHash()) == 1)
        {
            individual = std::move(UnitGenes());
        }
        mPopControl1.insert(individual.getHash());
        mPopulation1.push_back(std::move(individual));
        do
        {
            individual = std::move(UnitGenes());
        } while(mPopControl2.count(individual.getHash()) == 1);
        mPopControl2.insert(individual.getHash());
        mPopulation2.push_back(std::move(individual));
    }
}

template<class Race1, class Race2>
void UnitOptimizer<Race1,Race2>::mutate(size_t const N)
{

    auto mt = [=] (vector<UnitGenes>* population, unordered_set<double>* popControl, vector<UnitGenes>* newGenes)
    {
        size_t const nindividuals = std::min(N, population->size());
        std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<size_t> dist(0,nindividuals-1);

        auto chooseIndividual = std::bind(dist,gen);
        size_t number = static_cast<size_t>(std::round(mMutationRate*nindividuals));
        unordered_set<double> positions;

        do
        {
            positions.insert(chooseIndividual());
        } while(positions.size() < number);
        if(newGenes->capacity() < number)
        {
            newGenes->reserve(number);
        }
        for(size_t i : positions)
        {
            float const rate = 0.2;
            UnitGenes individual(population->at(i), rate);
            size_t count = 0;
            for(; count < 10 && popControl->count(individual.getHash()) == 1; ++count)
            {
                individual = std::move(UnitGenes(population->at(i), rate));
            }
            if(count >= 10)
            {
                do
                {
                    individual = std::move(UnitGenes());
                } while(popControl->count(individual.getHash()) == 1);
            }

            newGenes->push_back(std::move(individual));
        }

    };

    vector<UnitGenes> newGenes1, newGenes2;
    future<void> res = std::async(std::launch::async, mt, &mPopulation1, &mPopControl1, &newGenes1);
    //mt(&mPopulation1,&mPopControl1, &newGenes1);
    mt(&mPopulation2,&mPopControl2, &newGenes2);
    res.get();
    evaluate(newGenes1,mSim1);
    evaluate(newGenes2,mSim2);
    for(UnitGenes& genes : newGenes1)
    {
        mPopControl1.insert(genes.getHash());
        mPopulation1.push_back(std::move(genes));
    }
    for(UnitGenes& genes : newGenes2)
    {
        mPopControl2.insert(genes.getHash());
        mPopulation2.push_back(std::move(genes));
    }
}


template<class Race1, class Race2>
void UnitOptimizer<Race1,Race2>::crossover(size_t const N)
{
    auto co = [=] (vector<UnitGenes>* population, unordered_set<double>* popControl, vector<UnitGenes>* newGenes)
    {
        size_t const nindividuals = std::min(N, population->size());
        std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<size_t> dist(0,nindividuals-1);

        auto chooseIndividual = std::bind(dist,gen);
        size_t number = static_cast<size_t>(std::round(mReproductionRate*nindividuals));
        unordered_set<double> positions;
        if(newGenes->capacity() < number)
        {
            newGenes->reserve(number);
        }
        for(size_t i = 0; i < number; ++i)
        {
            do
            {
                positions.insert(chooseIndividual());
            } while(positions.size() < 4);

            long maxDist = 0, dist = 0;
            size_t pos1 = 0, pos2 = population->size()-1;
            for(size_t i : positions)
            {
                for(size_t j : positions)
                {
                    if(i != j)
                    {
                        dist = population->at(i).computeDistance(population->at(j));
                        if(dist > maxDist)
                        {
                            maxDist = dist;
                            pos1 = i;
                            pos2 = j;
                        }
                    }
                }
            }
            UnitGenes individual(population->at(pos1), population->at(pos2));
            size_t count = 0;
            for(; count < 10 && popControl->count(individual.getHash()) == 1; ++count)
            {
                individual = std::move(UnitGenes(population->at(pos1), population->at(pos2)));
            }
            if(count >= 10)
            {
                do
                {
                    individual = std::move(UnitGenes());
                } while(popControl->count(individual.getHash()) == 1);
            }
            newGenes->push_back(std::move(individual));
            positions.clear();
        }
    };
    vector<UnitGenes> newGenes1, newGenes2;
    future<void> res = std::async(std::launch::async, co, &mPopulation1, &mPopControl1, &newGenes1);
    //co(&mPopulation1, &mPopControl1, &newGenes1);
    co(&mPopulation2, &mPopControl2, &newGenes2);
    res.get();
    evaluate(newGenes1,mSim1);
    evaluate(newGenes2,mSim2);
    for(UnitGenes& genes : newGenes1)
    {
        mPopControl1.insert(genes.getHash());
        mPopulation1.push_back(std::move(genes));

    }
    for(UnitGenes& genes : newGenes2)
    {
        mPopControl2.insert(genes.getHash());
        mPopulation2.push_back(std::move(genes));
    }
}

template<class Race1, class Race2>
void UnitOptimizer<Race1,Race2>::resetPopulation()
{
    auto reset = [](UnitGenes& genes)
    {
        genes.fitness = Fitness();
        genes.count = 0;
    };

    std::for_each(mPopulation1.begin(), mPopulation1.end(), reset);
    std::for_each(mPopulation2.begin(), mPopulation2.end(), reset);
}

template<class Race1, class Race2>
void UnitOptimizer<Race1,Race2>::optimize(size_t iterations, size_t stepsPerIteration)
{

    mSurvival1.clear();
    mSurvival1.resize(iterations);
    mSurvival1.front() = 0.0;
    mSurvival2.clear();
    mSurvival2.resize(iterations);
    mSurvival2.front() = 0.0;
    auto select = [=] (vector<UnitGenes>* population, unordered_set<double>* popControl, float selectionRate)
    {
        std::sort(population->begin(),population->end());
        size_t const newSize = static_cast<size_t>(std::round(selectionRate*population->size()));
        while(population->size() > newSize)
        {
              popControl->erase(population->back().getHash());
              population->pop_back();
        }
    };


    std::cout << "Optimization started" << std::endl << "Number of forward simulations used: " << NTHREADS << std::endl;
    size_t minSize = 0;
    for(size_t i = 0; i < iterations; ++i)
    {
        mOldOptimum.first = mPopulation1.front();
        mOldOptimum.second = mPopulation2.front();
        std::cout << "Iteration " << i+1 << " started" << std::endl;
        initialize(mInitialPopulationSize);

        evaluate(mPopulation1,mSim1);
        evaluate(mPopulation2,mSim2);


        future<void> res = async(std::launch::async, select, &mPopulation1, &mPopControl1, mSelectionRate/2.);
        select(&mPopulation2, &mPopControl2, mSelectionRate/2.);
        res.get();

        for(size_t j = 0; j < stepsPerIteration; ++j)
        {
            pair<size_t,size_t> const popSize = std::make_pair<size_t, size_t>(mPopulation1.size(), mPopulation2.size());
            crossover(std::min(popSize.first, popSize.second));
            mutate(std::min(popSize.first, popSize.second));
            size_t const limit = std::pow(static_cast<size_t>(std::round(NTHREADS/mSelectionRate)),2);
            if(mPopulation1.size() < limit || mPopulation2.size() < limit)
            {
                break;
            }
            future<void> res = async(std::launch::async, select, &mPopulation1, &mPopControl1, mSelectionRate);
            select(&mPopulation2, &mPopControl2, mSelectionRate);
            res.get();
            std::cout << "Progress: " << static_cast<double>(j)/stepsPerIteration*100 << "%";
            std::cout << '\t' << "Population Size: " << mPopulation1.size() << '\r' << std::flush;
            printf("%c[2K", 27);

        }
        mOptFlag = true;
        minSize = std::min(mPopulation1.size(), mPopulation2.size());
        for(size_t i = 0; i < std::min(minSize, NTHREADS); ++i)
        {
            mSim1[i].setPlayer2Genes(mPopulation2[i]);
            mSim2[i].setPlayer2Genes(mPopulation1[i]);
        }

        std::cout << std::endl;

        printNBest(std::min(size_t(3),minSize));
	for(size_t k = 0; k < mPopulation1.size(); ++k)
	{
	    mSurvival1[i] += mPopulation1[k].count;
	}
	mSurvival1[i] /= mPopulation1.size();
	for(size_t k = 0; k < mPopulation2.size(); ++k)
	{
	    mSurvival2[i] += mPopulation2[k].count;
	}
	mSurvival2[i] /= mPopulation2.size();
        for(size_t i = 0; i < minSize; ++i)
        {
            ++mPopulation1[i].count;
            ++mPopulation2[i].count;
        }
//        long dist1 = mPopulation1.front().computeDistance(mOldOptimum.first), dist2 = mPopulation2.front().computeDistance(mOldOptimum.second);
//        if(dist1 < DELTA && dist2 < DELTA)
//        {
//            std::cout << "The distances between the optima of iteration " << i+1 << " and " << i << " are for both populations smaller than " << DELTA << " (" << dist1 << ", " << dist2 << ")" << std::endl;
//            break;
//        }

    }
    std::cout << "The individuals that survived the most iterations are:" << std::endl;
    printHighestCounts(10);
    std::cout << "Optimization finished" << std::endl;

}

template<class Race1, class Race2>
void UnitOptimizer<Race1,Race2>::optimize(size_t iterations, size_t stepsPerIteration, size_t initialPopulationSize, float selectionRate, float reproductionRate, float mutationRate)
{
    mInitialPopulationSize = initialPopulationSize;
    mSelectionRate = selectionRate;
    mReproductionRate = reproductionRate;
    mMutationRate = mutationRate;

    optimize(iterations, stepsPerIteration);
}

template<class Race1, class Race2>
void UnitOptimizer<Race1,Race2>::printBest()
{
    if(mPopulation1.empty() || mPopulation2.empty())
    {
        throw std::runtime_error("UnitOptimizer::printBest: The population is empty");
    }
    std::cout << "Player 1:" << std::endl << mPopulation1.front();
    std::cout << "Player 2:" << std::endl << mPopulation2.front();
}


template<class Race1, class Race2>
void UnitOptimizer<Race1,Race2>::printNBest(size_t const k)
{
    if(mPopulation1.empty() || mPopulation2.empty())
    {
        throw std::runtime_error("UnitOptimizer::printNBest: The population is empty");
    }
    size_t n = std::min(k,std::min(mPopulation1.size(),mPopulation2.size()));
    std::cout << "Player 1:" << std::endl;
    for(size_t i = 0; i < n; ++i)
    {
        std::cout << mPopulation1[i] << std::endl;
    }
    std::cout << "Player 2:" << std::endl;
    for(size_t i = 0; i < n; ++i)
    {
        std::cout << mPopulation2[i] << std::endl;
    }
}

template <class Race1, class Race2>
pair<UnitGenes,UnitGenes> UnitOptimizer<Race1,Race2>::getOptimum()
{
    if(mPopulation1.empty() || mPopulation2.empty())
    {
        throw std::runtime_error("UnitOptimizer::getOptimum: The population is empty");
    }
    return pair<UnitGenes,UnitGenes>(mPopulation1.front(),mPopulation2.front());
}

template <class Race1, class Race2>
pair<UnitGenes,UnitGenes> UnitOptimizer<Race1,Race2>::getOldOptimum()
{
    if(mOptFlag == false)
    {
        throw std::runtime_error("UnitOptimizer::getOldOptimum: optimize() has to be executed first");
    }
    return mOldOptimum;
}

template <class Race1, class Race2>
void UnitOptimizer<Race1,Race2>::printHighestCounts(size_t N)
{
    pair<vector<UnitGenes>,vector<UnitGenes>> res = std::move(getHighestCounts(N));

    size_t n = std::min(res.first.size(),res.second.size());
    std::cout << "Player 1:" << std::endl;
    for(size_t i = 0; i < n; ++i)
    {
        std::cout << res.first[i] << std::endl;
    }
    std::cout << "Player 2:" << std::endl;
    for(size_t i = 0; i < n; ++i)
    {
        std::cout << res.second[i] << std::endl;
    }
}

template <class Race1, class Race2>
pair<vector<UnitGenes>,vector<UnitGenes>> UnitOptimizer<Race1, Race2>::getHighestCounts(size_t N)
{
    vector<UnitGenes> pop1(mPopulation1);
    vector<UnitGenes> pop2(mPopulation2);
    auto cmp = [] (UnitGenes const& A, UnitGenes const& B)
    {
        return A.count > B.count;
    };

    std::sort(pop1.begin(), pop1.end(), cmp);
    std::sort(pop2.begin(), pop2.end(), cmp);
    pop1.resize(N);
    pop2.resize(N);
    return std::make_pair(pop1, pop2);
}
template <class Race1, class Race2>
pair<double,double> UnitOptimizer<Race1, Race2>::getHighscores() const
{
    return pair<double,double>(mPopulation1.front().fitness.score, mPopulation2.front().fitness.score);
}
template <class Race1, class Race2>
pair<vector<double>,vector<double>> UnitOptimizer<Race1, Race2>::getSurvivalAverage() const
{
    return pair<vector<double>,vector<double>>(mSurvival1, mSurvival2);
}

