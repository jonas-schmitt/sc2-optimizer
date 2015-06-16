#ifndef _SOGA_
#define _SOGA_

#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <algorithm>
#include <utility>
#include <functional>

#include "Chromosome.h"
#include "MicroSimulation.h"
#include "Utilities.h"

using std::vector;
using std::mt19937;
using std::bernoulli_distribution;
using std::string;
using std::pair;
using std::function;
using std::cout;
using std::endl;

struct Statistics
{
    double mean;
    double max;
    double sum;
    double var;
    Individual const* optimum;
    void print()
    {
        cout << "Total: " << sum << endl;
        cout << "Average: " << mean << endl;
        cout << "Maximum: " << max << endl;
        cout << "Variance: " << var << endl;
    }
};

template <typename T, typename U>
class SOGA final
{
private:
    double mutationProbability = 0.02;

    Statistics stats;

    size_t popSize;

    size_t NGenes;

    size_t NCrossoverPoints = 3;

    mt19937 generator;
    bernoulli_distribution dist1;
    std::uniform_int_distribution<size_t> dist2;
    std::uniform_real_distribution<double> dist3;
    std::uniform_int_distribution<size_t> dist4;
    bernoulli_distribution mutationDist;

    vector<Individual> pop;

    MicroSimulation<T,U> sim;


    vector<string> selectionFuncNames = {"Tournament Selection", "Roulette Wheel Selection", "Stochastic Universal Sampling"};
    vector<string> crossoverFuncNames = {"Single-Point Crossover", "Two-Point Crossover", "N-Point Crossover", "Uniform Crossover", "Three Parent Crossover"};
    vector<string> mutationFuncNames = {"Bit Flipping Mutation", "Interchanging Mutation", "Reversing Mutation"};

    size_t selectionChoice = 0;
    size_t crossoverChoice = 0;
    size_t mutationChoice = 0;



    function<vector<Individual *>(size_t const)> tournamentSelection = [&](size_t const N)
    {
        auto func = [&] ()
        {
            auto pickIndividual = std::bind(dist2, generator);
            size_t pos1 = pickIndividual();
            size_t pos2 = pickIndividual();
            Individual const& ind1 = pop[pos1];
            Individual const& ind2 = pop[pos2];
            return ind2.fitness.score < ind1.fitness.score ? pos1 : pos2;
        };

        vector<Individual *> res;
        res.reserve(N);
        for(size_t i = 0; i < N; ++i)
        {
            res.push_back(&pop[func()]);
        }
        return res;
    };


    function<vector<Individual *>(size_t const)> rouletteWheelSelection = [&](size_t const N)
    {
        auto func = [&] ()
        {
            auto spinWheel = std::bind(dist3, generator);
            auto cmp = [] (Individual const& ind, double val)
            {
                return ind.cdf < val;
            };

            auto low = std::lower_bound(pop.begin(), pop.end(), spinWheel(), cmp);
            return low;
        };

        vector<Individual *> res;
        res.reserve(N);
        for(size_t i = 0; i < N; ++i)
        {
            auto it = func();
            res.push_back(&(*it));
        }
        return res;
    };

    function<vector<Individual *>(size_t const)> stochasticUniversalSampling = [&](size_t const N)
    {
        auto func = [&] (double const p)
        {
            auto spinWheel = std::bind(dist3, generator);
            auto cmp = [] (Individual const& ind, double val)
            {
                return ind.cdf < val;
            };

            auto low = std::lower_bound(pop.begin(), pop.end(), p, cmp);
            return low;
        };
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
            auto it = func(p);
            res.push_back(&(*it));
        }
        return res;
    };

    // Crossover methods

    function<pair<Individual, Individual>(vector<Individual> const&)> singlePointCrossover = [&] (vector<Individual> const& parents)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
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
    };

    function<pair<Individual, Individual>(vector<Individual> const&)> twoPointCrossover = [&] (vector<Individual> const& parents)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
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
    };

    function<pair<Individual, Individual>(vector<Individual> const&)> nPointCrossover = [&] (vector<Individual> const& parents)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
        auto chooseBit = std::bind(dist4, generator);
        vector<size_t> bitPosArray(NCrossoverPoints);
        for(size_t& bitPos : bitPosArray)
        {
            bitPos = chooseBit();
        }
        std::sort(bitPosArray.begin(), bitPosArray.end());
        Individual child1(NGenes), child2(NGenes);
        auto swap = [] (Individual& A, Individual& B)
        {
            Individual& tmp = A;
            A = B;
            B = A;
        };
        Individual& A = child1;
        Individual& B = child2;

        size_t start = 0;
        for(size_t const bitPos : bitPosArray)
        {
            size_t const genePos = bitPos / NBITS;
            size_t const coPoint = bitPos - genePos * NBITS;
            for(size_t i = start; i < genePos; ++i)
            {
                A.chromosome[i] = parent1.chromosome[i];
                B.chromosome[i] = parent2.chromosome[i];
            }
            for(size_t i = 0; i < coPoint; ++i)
            {
                A.chromosome[genePos][i] = parent1.chromosome[genePos][i];
                B.chromosome[genePos][i] = parent2.chromosome[genePos][i];
            }
            swap(A,B);
            for(size_t i = coPoint; i < NBITS; ++i)
            {
                A.chromosome[genePos][i] = parent1.chromosome[genePos][i];
                B.chromosome[genePos][i] = parent2.chromosome[genePos][i];
            }
            start = genePos + 1;
        }
        for(size_t i = start; i < NGenes; ++i)
        {
            A.chromosome[i] = parent1.chromosome[i];
            B.chromosome[i] = parent2.chromosome[i];
        }
        return std::make_pair(child1, child2);
    };

    function<pair<Individual, Individual>(vector<Individual> const&)> uniformCrossover = [&] (vector<Individual> const& parents)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);

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
    };

    function<pair<Individual, Individual>(vector<Individual> const&)> threeParentCrossover = [&] (vector<Individual> const& parents)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
        Individual const& parent3 = parents.at(2);
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
        return std::make_pair(child, child);
    };

    function<void(Individual&)> bitFlipMutation = [&] (Individual& individual)
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
    };

    function<void(Individual&)> interchangingMutation = [&] (Individual& individual)
    {
        auto flipBiasedCoin = std::bind(mutationDist, generator);
        vector<pair<size_t, size_t>> positions;
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBITS; ++j)
            {
                if(flipBiasedCoin())
                {
                    positions.emplace_back(i,j);
                }
            }
        }
        std::shuffle(positions.begin(), positions.end(), generator);
        for(size_t i = 0; i < positions.size(); i += 2)
        {
            auto const& pos1 = positions[i];
            auto const& pos2 = positions[i+1];
            auto tmp = individual.chromosome[pos1.first][pos1.second];
            individual.chromosome[pos1.first][pos1.second] = individual.chromosome[pos2.first][pos2.second];
            individual.chromosome[pos2.first][pos2.second] = tmp;
        }
    };
    function<void(Individual&)> reversingMutation = [&] (Individual& individual)
    {
        auto flipCoin = std::bind(dist1, generator);
        auto flipBiasedCoin = std::bind(mutationDist, generator);
        vector<pair<size_t, size_t>> positions;
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBITS; ++j)
            {
                if(flipBiasedCoin())
                {
                    positions.emplace_back(i,j);
                }
            }
        }
        for(auto const& pos1 : positions)
        {
            pair<size_t, size_t> pos2;
            if((flipCoin() && !(pos1.first == NGenes-1 && pos1.second == NBITS-1)) || (pos1.first == 0 && pos1.second == 0))
            {
                // swap with right neighbour
                if(pos1.second == NBITS-1)
                {
                    pos2.first = pos1.first+1;
                    pos2.second = 0;
                }
                else
                {
                    pos2.first = pos1.first;
                    pos2.second = pos1.second + 1;
                }
            }
            else
            {
                if(pos1.second == 0)
                {
                    pos2.first = pos1.first-1;
                    pos2.second = NBITS-1;
                }
                else
                {
                    pos2.first = pos1.first;
                    pos2.second = pos1.second-1;
                }
            }
            auto tmp = individual.chromosome[pos1.first][pos1.second];
            individual.chromosome[pos1.first][pos1.second] = individual.chromosome[pos2.first][pos2.second];
            individual.chromosome[pos2.first][pos2.second] = tmp;
        }
    };


    vector<function<vector<Individual *>(size_t const) > > selectionFuncs = {tournamentSelection, rouletteWheelSelection, stochasticUniversalSampling};
    vector<function<pair<Individual, Individual>(vector<Individual> const&)> > crossoverFuncs = {singlePointCrossover, twoPointCrossover, nPointCrossover, uniformCrossover, threeParentCrossover};
    vector<function<void(Individual &) > > mutationFuncs = {bitFlipMutation, interchangingMutation, reversingMutation};





    void evaluate(Individual& ind)
    {
        sim.setPlayer1Chromosome(ind.chromosome);
        ind.fitness = std::move(sim.run(true));
    }

    void setGoal(Chromosome const& goal)
    {
        sim.setPlayer2Chromosome(goal);
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
                stats.optimum = &ind;
            }
            stats.sum += ind.fitness.score;
        }
        stats.mean = stats.sum / pop.size();

        stats.var = 0.0;
        for(Individual const & ind : pop)
            stats.var += (stats.mean-ind.fitness.score)*(stats.mean-ind.fitness.score);
        stats.var /= pop.size();

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

        pop.reserve(2*popSize);
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
                stats.optimum = &pop[i];
            }
            stats.sum += pop[i].fitness.score;
        }
        stats.mean = stats.sum / popSize;
        computeCDF();

    }

    void optimize(Chromosome const& goal, size_t const iterations)
    {
        mutationDist = std::move(std::bernoulli_distribution(mutationProbability));
        setGoal(goal);
        computeStatistics();
        computeCDF();
        bool const tpco = crossoverChoice == crossoverFuncs.size()-1;
        for(size_t i = 0; i < iterations; ++i)
        {
            // apply genetic algorithm

            vector<Individual *> selected;
            size_t add = 2;
            if(tpco)
            {
                selected = selectionFuncs[selectionChoice](3*popSize);
                add = 3;
            }
            else
            {
                selected = selectionFuncs[selectionChoice](popSize);
            }
            for(size_t i = 0; i < selected.size(); i += add)
            {
                pair<Individual, Individual> children;
                if(tpco)
                {
                    children = crossoverFuncs[crossoverChoice]({*selected[i], *selected[i+1], *selected[i+2]});
                    mutationFuncs[mutationChoice](children.first);
                    evaluate(children.first);
                    pop.push_back(children.first);
                }
                else
                {
                    children = crossoverFuncs[crossoverChoice]({*selected[i], *selected[i+1]});
                    mutationFuncs[mutationChoice](children.first);
                    evaluate(children.first);

                    mutationFuncs[mutationChoice](children.second);
                    evaluate(children.second);

                    pop.push_back(children.first);
                    pop.push_back(children.second);
                }
            }
            auto cmp = [] (Individual const& ind1, Individual const& ind2)
            {
                return ind1.fitness.score > ind2.fitness.score;
            };
            sort(pop.begin(), pop.end(), cmp);
            pop.resize(popSize);

            computeStatistics();
            computeCDF();
        }

    }

    size_t getNumberOfSelectionOperators()
    {
        return selectionFuncs.size();
    }

    size_t getNumberOfCrossoverOperators()
    {
        return crossoverFuncs.size();
    }

    size_t getNumberOfMutationOperators()
    {
        return mutationFuncs.size();
    }

    Statistics getStatistics()
    {
        return stats;
    }

    void setSelection(size_t const value)
    {
        selectionChoice = value < selectionFuncs.size() ? value : selectionChoice;
    }

    void setCrossover(size_t const value)
    {
        crossoverChoice = value < crossoverFuncs.size() ? value : crossoverChoice;
    }

    void setMutation(size_t const value)
    {
        mutationChoice = value < mutationFuncs.size() ? value : mutationChoice;
    }

    string getSelectionOperatorName()
    {
        return selectionFuncNames[selectionChoice];
    }
    string getCrossoverOperatorName()
    {
        return crossoverFuncNames[selectionChoice];
    }
    string getMutationOperatorName()
    {
        return mutationFuncNames[selectionChoice];
    }




};


#endif // SOGA

