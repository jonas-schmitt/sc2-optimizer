#ifndef _SOGA_
#define _SOGA_

#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <algorithm>
#include <utility>

#include "Chromosome.h"
#include "MicroSimulation.h"
#include "Utilities.h"

using std::vector;
using std::mt19937;
using std::bernoulli_distribution;
using std::string;
using std::pair;

struct Statistics
{
    double avg;
    double max;
    double sum;
    double var;
};

template <typename T, typename U>
class SOGA final
{
private:
    double mutationProbability = 0.02;

    Statistics stats;

    size_t popSize;

    size_t NGenes;

    mt19937 generator;
    bernoulli_distribution dist1;
    std::uniform_int_distribution<size_t> dist2;
    std::uniform_real_distribution<double> dist3;
    std::uniform_int_distribution<size_t> dist4;
    bernoulli_distribution mutationDist;

    vector<Individual> pop;

    MicroSimulation<T,U> sim;

    Individual const* optimum;


    // Selection Methods
    Individual& tournamentSelect()
    {
        auto pickIndividual = std::bind(dist2, generator);
        Individual& ind1 = pop[pickIndividual()];
        Individual& ind2 = pop[pickIndividual()];
        return ind2.fitness.score < ind1.fitness.score ? ind1 : ind2;
    }


    Individual& rouletteWheelSelection(double const p)
    {
        auto cmp = [] (Individual const& ind, double val)
        {
            return ind.cdf < val;
        };

        auto low = std::lower_bound(pop.begin(), pop.end(), p);
        return *low;
    }

    Individual& rouletteWheelSelection()
    {
        auto spinWheel = std::bind(dist3, generator);
        return rouletteWheelSelection(spinWheel());
    }

    vector<Individual *> tournamentSelection(size_t const N)
    {
        vector<Individual *> res;
        res.reserve(N);
        for(size_t i = 0; i < N; ++i)
        {
            res.push_back(&tournamentSelection());
        }
        return res;
    }


    vector<Individual *> rouletteWheelSelection(size_t const N)
    {
        vector<Individual *> res;
        res.reserve(N);
        for(size_t i = 0; i < N; ++i)
        {
            res.push_back(&rouletteWheelSelection());
        }
        return res;
    }

    vector<Individual *> stochasticUniversalSampling(size_t const N)
    {
        vector<Individual *> res;
        res.reserve(N);
        double const distance = stats.sum / N;
        std::uniform_real_distribution<double> distribution(0,distance);
        double const start = distribution(generator);
        vector<double> pointers;
        pointers.reserve(N);
        for(size_t i = 0; i < N; ++i)
        {
            pointers.push_back(start + i * distance);
        }
        for(double const p : pointers)
        {
            res.push_back(&rouletteWheelSelection(p));
        }
        return res;
    }

    // Crossover methods

    pair<Individual, Individual> singlePointCrossover(Individual const& parent1, Individual const &parent2)
    {
        auto chooseBit = std::bind(dist4, generator);
        size_t bitPos = chooseBit();
        size_t const genePos = bitPos / NBITS;
        size_t const coPoint = bitPos - genePos * NBITS;
        Individual child1(NGenes), child2(NGenes);
        for(size_t i = 0; i < genePos; ++i)
        {
            child1.chromosome[i] = parent1.chromosome[i];
            child2.chromosome[i] = parent2.chromosome[i];
        }
        for(size_t i = 0; i < coPoint; ++i)
        {
            child1.chromosome[genePos][i] = parent1.chromosome[genePos][i];
            child2.chromosome[genePos][i] = parent2.chromosome[genePos][i];
        }
        for(size_t i = coPoint; i < NBITS; ++i)
        {
            child1.chromosome[genePos][i] = parent2.chromosome[genePos][i];
            child2.chromosome[genePos][i] = parent1.chromosome[genePos][i];
        }
        for(size_t i = genePos + 1; i < NGenes; ++i)
        {
            child1.chromosome[i] = parent2.chromosome[i];
            child2.chromosome[i] = parent1.chromosome[i];
        }
        return std::make_pair(child1, child2);
    }

    pair<Individual, Individual> twoPointCrossover(Individual const& parent1, Individual const &parent2)
    {
        auto chooseBit = std::bind(dist4, generator);
        size_t bitPos1 = chooseBit();
        size_t bitPos2 = chooseBit();
        while(bitPos1 == bitPos2)
        {
            bitPos2 = chooseBit();
        }

        if(bitPos1 > bitPos2)
        {
            size_t tmp = bitPos2;
            bitPos2 = bitPos1;
            bitPos1 = tmp;
        }
        size_t genePos1 = bitPos1 / NBITS;
        size_t coPoint1 = bitPos1 - genePos1 * NBITS;

        size_t genePos2 = bitPos2 / NBITS;
        size_t coPoint2 = bitPos2 - genePos2 * NBITS;

        Individual child1(NGenes), child2(NGenes);
        for(size_t i = 0; i < genePos1; ++i)
        {
            child1.chromosome[i] = parent1.chromosome[i];
            child2.chromosome[i] = parent2.chromosome[i];
        }
        for(size_t i = 0; i < coPoint1; ++i)
        {
            child1.chromosome[genePos1][i] = parent1.chromosome[genePos1][i];
            child2.chromosome[genePos1][i] = parent2.chromosome[genePos1][i];
        }
        for(size_t i = coPoint1; i < NBITS; ++i)
        {
            child1.chromosome[genePos1][i] = parent2.chromosome[genePos1][i];
            child2.chromosome[genePos1][i] = parent1.chromosome[genePos1][i];
        }
        for(size_t i = genePos1 + 1; i < genePos2; ++i)
        {
            child1.chromosome[i] = parent2.chromosome[i];
            child2.chromosome[i] = parent1.chromosome[i];
        }
        for(size_t i = 0; i < coPoint2; ++i)
        {
            child1.chromosome[genePos2][i] = parent2.chromosome[genePos2][i];
            child2.chromosome[genePos2][i] = parent1.chromosome[genePos2][i];
        }
        for(size_t i = coPoint2; i < NBITS; ++i)
        {
            child1.chromosome[genePos2][i] = parent1.chromosome[genePos2][i];
            child2.chromosome[genePos2][i] = parent2.chromosome[genePos2][i];
        }
        for(size_t i = genePos2 + 1; i < NGenes; ++i)
        {
            child1.chromosome[i] = parent1.chromosome[i];
            child2.chromosome[i] = parent2.chromosome[i];
        }

        return std::make_pair(child1, child2);
    }

    pair<Individual, Individual> uniformCrossover(Individual const& parent1, Individual const &parent2)
    {
        auto flipCoin = std::bind(dist1, generator);
        Individual child1(NGenes), child2(NGenes);
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBITS; ++j)
            {
                if(flipCoin())
                {
                    child1.chromosome[i][j] = parent1.chromosome[i][j];
                    child2.chromosome[i][j] = parent2.chromosome[i][j];
                }
                else
                {
                    child1.chromosome[i][j] = parent2.chromosome[i][j];
                    child2.chromosome[i][j] = parent1.chromosome[i][j];
                }
            }
        }
        return std::make_pair(child1, child2);
    }

    Individual threeParentCrossover(Individual const& parent1, Individual const &parent2, Individual const& parent3)
    {
        Individual child(NGenes);
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBITS; ++j)
            {

                if(parent1.chromosome[i][j] == parent2.chromosome[i][j])
                {
                    child.chromosome[i][j] = parent1.chromosome[i][j];
                }
                else
                {
                    child.chromosome[i][j] = parent3.chromosome[i][j];
                }

            }
        }
        return child;
    }

    void mutation(Individual& individual)
    {
        auto flipBiasedCoin = std::bind(mutationDist, generator);
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBITS; ++j)
            {
                if(flipBiasedCoin())
                {
                    individual.chromosome[i].flip(j);
                }
            }
        }
    }





    void evaluate(Individual& ind)
    {
        sim.setPlayer1Chromosome(ind.chromosome);
        ind.fitness = std::move(sim.run(true));
    }

    void setGoal(Chromosome const& goal)
    {
        sim.setPlayer2Chromosome(goal.size());
    }

    void computeCDF()
    {
        double const inv_sum = 1.0/stats.sum;
        double tmp = 0.0;
        for(Individual& ind : pop)
        {
            tmp += ind.fitness.score * inv_sum;
            ind.cdf = tmp;
        }
    }

    void computeStatistics()
    {
        stats.max = 0.0;
        stats.sum = 0.0;
        for(Individual& ind : pop)
        {
            if(ind.fitness.score > stats.max)
            {
                stats.max = ind.fitness.score;
                optimum = &ind;
            }
            stats.sum += ind.fitness.score;
        }
        stats.avg = stats.sum / pop.size();
    }


public:

    SOGA(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2, size_t popSize, vector<string> const & buildList1, vector<string> const & buildList2)
        :  popSize(popSize), generator(std::chrono::system_clock::now().time_since_epoch().count()), dist1(0.5), dist2(0,popSize-1), dist3(0,1.0),
          sim(minPos, maxPos, filePath1, filePath2)
    {
        sim.initBothPlayers(buildList1, buildList2);
        NGenes = sim.getPlayer1ChromosomeLength();
        dist4 = std::uniform_int_distribution<size_t>(1, NGenes*NBITS - 1);
        Chromosome initChrom(NGenes);
        for(auto& gene : initChrom)
        {
            gene.flip(gene.size()-1);
        }

        sim.setPlayer2Chromosome(initChrom);

        pop.clear();
        pop.resize(popSize);


        stats.sum = 0.0;
        stats.max = 0.0;

        auto flipCoin = std::bind(dist1, generator);
        for(size_t i = 0; i < popSize; ++i)
        {
            Chromosome& chrom = pop[i].chromosome;
            chrom.resize(NGenes);
            for(auto& gene : chrom)
            {
                for(size_t j = 0; j < gene.size(); ++j)
                {
                    gene.set(j, flipCoin());
                }
            }
            evaluate(pop[i]);
            if(pop[i].fitness.score > stats.max)
            {
                stats.max = pop[i].fitness.score;
                optimum = &pop[i];
            }
            stats.sum += pop[i].fitness.score;
        }
        stats.avg = stats.sum / popSize;
        computeCDF();

    }

    void optimize(Chromosome const& goal, size_t const iterations)
    {
        mutationDist = std::move(std::bernoulli_distribution(mutationProbability));
        setGoal(goal);
        for(size_t i = 0; i < iterations; ++i)
        {
            // TODO apply genetic algorithm
            computeStatistics();
            computeCDF();
        }

    }


};


#endif // SOGA

