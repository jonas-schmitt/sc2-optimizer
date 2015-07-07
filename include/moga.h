#ifndef _MOGA_
#define _MOGA_

#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <algorithm>
#include <utility>
#include <functional>
#include <unordered_set>
#include <thread>
#include <unistd.h>
#include <omp.h>

#include "Chromosome.h"
#include "MicroSimulation.h"
#include "Utilities.h"

using std::vector;
using std::mt19937;
using std::bernoulli_distribution;
using std::string;
using std::pair;
using std::unordered_set;
using std::function;
using std::cout;
using std::endl;


template <typename T, typename U, Player const player>
class MOGA final
{
private:
    double mutationProbability = 0.01;

    Statistics stats;

    size_t popSize;

    size_t NGenes;

    size_t NCrossoverPoints;

    mt19937 generator;
    bernoulli_distribution flipCoin;
    std::uniform_int_distribution<size_t> chooseIndividual;
    std::uniform_real_distribution<double> spinWheel;
    std::uniform_int_distribution<size_t> chooseBit;
    bernoulli_distribution mutationDist;

    vector<Individual> pop;

    vector<vector<MicroSimulation<T,U>>> sim;

    unordered_set<size_t> control;

    vector<string> selectionFuncNames = {"Tournament Selection", "Roulette Wheel Selection", "Stochastic Universal Sampling"};
    vector<string> crossoverFuncNames = {"Single-Point Crossover", "Two-Point Crossover", "N-Point Crossover", "Uniform Crossover"};
    vector<string> mutationFuncNames = {"Bit Flipping Mutation", "Interchanging Mutation", "Reversing Mutation"};

    size_t selectionChoice = 0;
    size_t crossoverChoice = 0;
    size_t mutationChoice = 0;

    double onlinePerformance = 0.0;
    double offlinePerformance = 0.0;

    size_t NGoals;



    function<vector<Individual *>(size_t const, mt19937&, vector<Individual>&)> tournamentSelection = [&](size_t const N, mt19937& generator, vector<Individual>& pop)
    {
        std::uniform_int_distribution<size_t> dist(0,pop.size()-1);
        auto func = [&] (std::uniform_int_distribution<size_t>& dist, mt19937& generator)
        {
            size_t const pos1 = dist(generator);
            size_t pos2 = dist(generator);
            while(pos2 == pos1) pos2 = dist(generator);
            Individual const& lhs = pop[pos1];
            Individual const& rhs = pop[pos2];
            if(lhs.rank < rhs.rank)
            {
                return pos1;
            }
            else if(lhs.rank == rhs.rank && lhs.distance > rhs.distance)
            {
                return pos1;
            }
            else
            {
                return pos2;
            }
        };

        vector<Individual *> res;
        res.reserve(N);
        for(size_t i = 0; i < N; ++i)
        {
            res.push_back(&pop[func(dist, generator)]);
        }
        return res;
    };


    function<vector<Individual *>(size_t const, mt19937&, vector<Individual>&)> rouletteWheelSelection = [&](size_t const N, mt19937& generator, vector<Individual>& pop)
    {
        std::uniform_real_distribution<double> dist(0,1.0);
        auto func = [&] (double const p)
        {
            auto cmp = [] (Individual const& ind, double val)
            {
                return ind.cdf < val;
            };

            auto low = std::lower_bound(pop.begin(), pop.end(), p, cmp);
            return low;
        };

        vector<Individual *> res;
        res.reserve(N);
        for(size_t i = 0; i < N; ++i)
        {
            auto it = func(dist(generator));
            res.push_back(&(*it));
        }
        return res;
    };

    function<vector<Individual *>(size_t const, mt19937&, vector<Individual>&)> stochasticUniversalSampling = [&](size_t const N, mt19937& generator, vector<Individual>& pop)
    {
        auto func = [&] (double const p)
        {
            auto cmp = [] (Individual const& ind, double val)
            {
                return ind.cdf < val;
            };

            auto low = std::lower_bound(pop.begin(), pop.end(), p, cmp);
            return low;
        };
        vector<Individual *> res;
        res.reserve(N);
        double const total = pop.front ().total;
        double const distance = total / N;
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
            auto it = func(p/total);
            res.push_back(&(*it));
        }
        std::shuffle(res.begin(), res.end(), generator);
        return res;
    };

    // Crossover methods

    function<pair<Individual, Individual>(vector<Individual> const&, mt19937&, size_t const, size_t const)> singlePointCrossover = [&] (vector<Individual> const& parents, mt19937& generator, size_t const NGenes, size_t const NBits)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
        std::uniform_int_distribution<size_t> dist(1,NGenes*NBits-2);
        size_t bitPos = dist(generator);
        size_t const genePos = bitPos / NBits;
        size_t const coPoint = bitPos - genePos * NBits;
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

    function<pair<Individual, Individual>(vector<Individual> const&, mt19937&, size_t const, size_t const)> twoPointCrossover = [&] (vector<Individual> const& parents, mt19937& generator, size_t const NGenes, size_t const NBits)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
        std::uniform_int_distribution<size_t> dist(1,NGenes*NBits-2);
        size_t bitPos1 = dist(generator);
        size_t bitPos2 = dist(generator);
        while(bitPos1 / NBits == bitPos2 / NBits)
        {
            bitPos2 = dist(generator);
        }

        if(bitPos1 > bitPos2)
        {
            size_t const tmp = bitPos2;
            bitPos2 = bitPos1;
            bitPos1 = tmp;
        }
        size_t const genePos1 = bitPos1 / NBits;
        size_t const coPoint1 = bitPos1 - genePos1 * NBits;

        size_t const genePos2 = bitPos2 / NBits;
        size_t const coPoint2 = bitPos2 - genePos2 * NBits;

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
        for(size_t i = coPoint1; i < NBits; ++i)
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
        for(size_t i = coPoint2; i < NBits; ++i)
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

    function<pair<Individual, Individual>(vector<Individual> const&, mt19937&, size_t const, size_t const)> nPointCrossover = [&] (vector<Individual> const& parents, mt19937& generator, size_t const NGenes, size_t const NBits)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
        std::uniform_int_distribution<size_t> dist(1,NGenes*NBits-2);

        vector<size_t> bitPosArray((NGenes -1) / 2);
        unordered_set<size_t> positions;
        for(size_t& bitPos : bitPosArray)
        {
            do {
            bitPos = dist(generator);
            } while(positions.count(bitPos / NBits) == 1);
            positions.insert(bitPos / NBits);
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
            size_t const genePos = bitPos / NBits;
            size_t const coPoint = bitPos - genePos * NBits;
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
            for(size_t i = coPoint; i < NBits; ++i)
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

    function<pair<Individual, Individual>(vector<Individual> const&, mt19937&, size_t const, size_t const)> uniformCrossover = [&] (vector<Individual> const& parents, mt19937& generator, size_t const NGenes, size_t const NBits)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);

        std::bernoulli_distribution dist(0.5);
        Individual child1(NGenes), child2(NGenes);
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBits; ++j)
            {
                if(dist(generator))
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

    function<pair<Individual, Individual>(vector<Individual> const&, mt19937&, size_t const, size_t const)> threeParentCrossover = [&] (vector<Individual> const& parents, mt19937& generator, size_t const NGenes, size_t const NBits)
    {
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
        Individual const& parent3 = parents.at(2);
        Individual child(NGenes);
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBits; ++j)
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

    function<void(Individual&, mt19937&, std::bernoulli_distribution&, size_t const, size_t const)> bitFlipMutation = [&] (Individual& individual, mt19937& generator, std::bernoulli_distribution& mutationDist, size_t const NGenes, size_t const NBits)
    {
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBits; ++j)
            {
                if(mutationDist(generator))
                {
                    individual.chromosome[i].flip(j);
                }
            }
        }
    };

    function<void(Individual&, mt19937&, std::bernoulli_distribution&, size_t const, size_t const)> interchangingMutation = [&] (Individual& individual, mt19937& generator, std::bernoulli_distribution& mutationDist, size_t const NGenes, size_t const NBits)
    {
        vector<pair<size_t, size_t>> positions;
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBits; ++j)
            {
                if(mutationDist(generator))
                {
                    positions.emplace_back(i,j);
                }
            }
        }
        if(positions.size() == 0) return;
        std::shuffle(positions.begin(), positions.end(), generator);
        for(size_t i = 0; i < positions.size()-1; i += 2)
        {
            auto const& pos1 = positions[i];
            auto const& pos2 = positions[i+1];
            auto tmp = individual.chromosome[pos1.first][pos1.second];
            individual.chromosome[pos1.first][pos1.second] = individual.chromosome[pos2.first][pos2.second];
            individual.chromosome[pos2.first][pos2.second] = tmp;
        }
    };
    function<void(Individual&, mt19937&, std::bernoulli_distribution&, size_t const, size_t const)> reversingMutation = [&] (Individual& individual, mt19937& generator, std::bernoulli_distribution& mutationDist, size_t const NGenes, size_t const NBits)
    {
        std::bernoulli_distribution coin(0.5);
        vector<pair<size_t, size_t>> positions;
        for(size_t i = 0; i < NGenes; ++i)
        {
            for(size_t j = 0; j < NBits; ++j)
            {
                if(mutationDist(generator))
                {
                    positions.emplace_back(i,j);
                }
            }
        }
        for(auto const& pos1 : positions)
        {
            pair<size_t, size_t> pos2;
            if((coin(generator) && !(pos1.first == NGenes-1 && pos1.second == NBits-1)) || (pos1.first == 0 && pos1.second == 0))
            {
                // swap with right neighbour
                if(pos1.second == NBits-1)
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
                    pos2.second = NBits-1;
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


    vector<function<vector<Individual *>(size_t const, mt19937&, vector<Individual>&) > > selectionFuncs = {tournamentSelection, rouletteWheelSelection, stochasticUniversalSampling};
    vector<function<pair<Individual, Individual>(vector<Individual> const&, mt19937&, size_t const, size_t const)> > crossoverFuncs = {singlePointCrossover, twoPointCrossover, nPointCrossover, uniformCrossover};
    vector<function<void(Individual&, mt19937&, std::bernoulli_distribution&, size_t const, size_t const) > > mutationFuncs = {bitFlipMutation, interchangingMutation, reversingMutation};





    void evaluate(vector<Individual>& pop)
    {

        #pragma omp parallel for schedule(runtime)
        for(size_t i = 0; i < pop.size(); ++i)
        {
            double const value = 1.0/NGoals;
            pop[i].fitness = 0;
            for(size_t j = 0; j < NGoals; ++j)
            {
                if(player == Player::first)
                {
                    sim[omp_get_thread_num()][j].setPlayer1Chromosome(pop[i].chromosome);
                    pop[i].fitness += sim[omp_get_thread_num()][j].run(true, Player::first);
                }
                else
                {
                    sim[omp_get_thread_num()][j].setPlayer2Chromosome(pop[i].chromosome);
                    pop[i].fitness += sim[omp_get_thread_num()][j].run(true, Player::second);
                }

            }
            pop[i].fitness *= value;
        }
    }

    void setGoals(vector<Chromosome> const& goals)
    {
        if(goals.size() != NGoals)
        {
            throw std::invalid_argument("MOGA::setGoals(): Invalid number of arguments");
        }
        #pragma omp parallel
        {
            for(size_t j = 0; j < NGoals; ++j)
            {
                if(player == Player::first)
                {
                    sim[omp_get_thread_num()][j].setPlayer2Chromosome(goals[j]);
                }
                else
                {
                    sim[omp_get_thread_num()][j].setPlayer1Chromosome(goals[j]);
                }
            }
            #pragma omp barrier
        }
    }

    void computeCDF()
    {
        double const inv_sum = 1.0/stats.sum;
        double tmp = 0.0;
        for(Individual& ind : pop)
        {
            tmp += ind.fitness.score * inv_sum;
            ind.cdf = tmp;
            ind.total = stats.sum;
        }
    }

    void computeStatistics()
    {
        vector<double> v(pop.size());
        for(size_t i = 0; i < v.size(); ++i)
        {
            v[i] = pop[i].fitness.score;
        }
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        double mean = sum / v.size();
        std::vector<double> diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(),
                       std::bind2nd(std::minus<double>(), mean));
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / v.size());
        stats.sum = sum;
        stats.mean = mean;
        stats.max = *std::max_element(v.begin(), v.end());
        stats.stdev = stdev;

    }

    void nondominatedSort()
    {
        // nondominated fronts
        vector<vector<size_t>> fronts(1);

        // determine the domination count and set of dominated individuals for all members of the population
        for(size_t i = 0; i < pop.size(); ++i)
        {
            pop[i].dominationCount = 0;
            pop[i].dominationSet.clear();
            pop[i].distance = 0.0;
            for(size_t j = 0; j < pop.size(); ++j)
            {
                if(i == j) continue;
                if(pop[i].dominates(pop[j]))
                {
                    pop[i].dominationSet.push_back(j);
                }
                else if(pop[j].dominates(pop[i]))
                {
                    ++pop[i].dominationCount;
                }
            }
            // if the individual is not dominated by any other individuals, include it in the first nondominated front
            if(pop[i].dominationCount == 0)
            {
                fronts[0].push_back(i);
                pop[i].rank = 1;
            }
        }

        // determine the individuals in each subsequent front
        size_t r = 1;
        size_t count = fronts[0].size();
        fronts.emplace_back();
        for(size_t i = 0; count < popSize && !fronts[i].empty(); ++i)
        {
            for(size_t const p : fronts[i])
            {
                for(size_t const q : pop[p].dominationSet)
                {
                    if(pop[q].dominationCount > 0)
                    {
                        if(--pop[q].dominationCount == 0)
                        {
                            fronts[i+1].push_back(q);
                            pop[q].rank = r + 1;
                        }
                    }
                }
            }
            count += fronts[i+1].size();
            fronts.emplace_back();
            ++r;
        }


        vector<Individual> newPop;
        newPop.reserve(popSize);
        fronts.pop_back();


        auto cmp_damage = [&] (size_t const p, size_t const q)
        {
            return pop[p].fitness.damage < pop[q].fitness.damage;
        };
        auto cmp_health = [&] (size_t const p, size_t const q)
        {
            return pop[p].fitness.health < pop[q].fitness.health;
        };

        auto cmp_time = [&] (size_t const p, size_t const q)
        {
            return pop[p].fitness.timeSteps < pop[q].fitness.timeSteps;
        };


        // diversity preservation
        for(vector<size_t>& front : fronts)
        {

            // compute crowding distance for damage
            std::sort(front.begin(), front.end(), cmp_damage);
            pop[front.front()].distance = INF;
            pop[front.back()].distance = INF;
            if(front.size() > 2)
            {
                for(size_t i = 1; i < front.size()-2; ++i)
                {
                    pop[front[i]].distance += pop[front[i+1]].fitness.damage - pop[front[i-1]].fitness.damage;
                }
            }



            // compute crowding distance for health
            std::sort(front.begin(), front.end(), cmp_health);
            pop[front.front()].distance = INF;
            pop[front.back()].distance = INF;
            if(front.size() > 2)
            {
                for(size_t i = 1; i < front.size()-2; ++i)
                {
                    if(pop[front[i]].distance < INF)
                    {
                        pop[front[i]].distance += pop[front[i+1]].fitness.health - pop[front[i-1]].fitness.health;
                    }
                }
            }

            // compute crowding distance for health
            std::sort(front.begin(), front.end(), cmp_time);
            pop[front.front()].distance = INF;
            pop[front.back()].distance = INF;
            if(front.size() > 2)
            {
                for(size_t i = 1; i < front.size()-2; ++i)
                {
                    if(pop[front[i]].distance < INF)
                    {
                        pop[front[i]].distance += pop[front[i+1]].fitness.timeSteps - pop[front[i-1]].fitness.timeSteps;
                    }
                }
            }


            // if the front does not fit completely in the population, sort the individuals according to the crowded comparison operator
            if(newPop.size() + front.size() > popSize)
            {
                auto cmp_distance = [&] (size_t const p, size_t const q)
                {
                    return pop[p].distance > pop[q].distance;
                };
                sort(front.begin(), front.end(), cmp_distance);
            }

            // consecutively include individuals from each front into the new population until it is filled
            for(size_t const p : front)
            {
                if(newPop.size() == popSize) break;
                newPop.push_back(pop[p]);
            }

        }
        // replace the old population with the new one
        pop = newPop;
    }


public:

    typedef T race1;
    typedef U race2;

    MOGA(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2, size_t popSize, vector<string> const & buildList1, vector<string> const & buildList2, size_t const nGoals)
        :  popSize(popSize), generator(std::chrono::system_clock::now().time_since_epoch().count()), flipCoin(0.5), chooseIndividual(0,popSize-1), spinWheel(0,1.0), NGoals(nGoals)
    {
//        sim.resize(omp_get_max_threads());
//        for(int i = 0; i < omp_get_max_threads(); ++i)
//        {
//            sim[i].reserve(NGoals);
//            for(size_t j = 0; j < NGoals; ++j)
//            {
//                sim[i].emplace_back(minPos, maxPos, filePath1, filePath2);
//                sim[i].back().initBothPlayers(buildList1, buildList2);
//            }

//        }

        sim.reserve(omp_get_max_threads());
        #pragma omp parallel
        {
            for(size_t i = 0; i < omp_get_max_threads(); ++i)
            {
                #pragma omp critical
                {
                    if(i == omp_get_thread_num())
                    {
                        sim.emplace_back();
                        sim[omp_get_thread_num()].reserve(NGoals);
                        for(size_t j = 0; j < NGoals; ++j)
                        {
                            sim[omp_get_thread_num()].emplace_back(minPos, maxPos, filePath1, filePath2);
                            sim[omp_get_thread_num()].back().initBothPlayers(buildList1, buildList2);
                        }
                    }
                }
                #pragma omp barrier
            }
        }

        if(player == Player::first)
        {
            NGenes = sim[0][0].getPlayer1ChromosomeLength();
            NCrossoverPoints = buildList1.size();
        }
        else
        {
            NGenes = sim[0][0].getPlayer2ChromosomeLength();
            NCrossoverPoints = buildList2.size();
        }
        chooseBit = std::uniform_int_distribution<size_t>(1, NGenes*NBITS - 1);
        vector<Chromosome> initChroms(NGoals);
        for(auto& chrom : initChroms)
        {
            chrom.resize(NGenes);
            for(auto& gene : chrom)
            {
                for(size_t i = 0; i < gene.size(); ++i)
                {
                    gene.set(i, flipCoin(generator));
                }
            }
        }
        setGoals(initChroms);


        pop.reserve(2*popSize);
        pop.resize(popSize);


        stats.sum = 0.0;
        stats.max = 0.0;

        for(size_t i = 0; i < popSize; ++i)
        {
            Chromosome& chrom = pop[i].chromosome;
            chrom.resize(NGenes);
            double hash;
            do
            {
                for(auto& gene : chrom)
                {
                    for(size_t j = 0; j < gene.size(); ++j)
                    {
                        gene.set(j, flipCoin(generator));
                    }
                }
                hash = pop[i].computeHash();
            }
            while(control.count(hash) == 1);
            control.insert(hash);
        }
        evaluate(pop);
        computeStatistics();
        computeCDF();

    }

    void optimize(vector<Chromosome> const& goals, size_t const iterations)
    {
        onlinePerformance = 0.0;
        offlinePerformance = 0.0;

        mutationDist = std::move(std::bernoulli_distribution(mutationProbability));
        setGoals(goals);
        evaluate(pop);

        computeStatistics();
        computeCDF();
        for(size_t i = 0; i < iterations; ++i)
        {
            // apply genetic algorithm

            vector<Individual> newPop;
            newPop.reserve(popSize);
            vector<Individual *> selected;

            selected = selectionFuncs[selectionChoice](popSize, generator, pop);

            for(size_t count = 0; count < 100; ++count)
            {
                if(selected.size() == 0) break; // Just for safety
                for(size_t i = 0; i < selected.size()-1 && newPop.size() < popSize; i += 2)
                {
                    pair<Individual, Individual> children;

                    children = crossoverFuncs[crossoverChoice]({*selected[i], *selected[i+1]}, generator, NGenes, NBITS);
                    mutationFuncs[mutationChoice](children.first, generator, mutationDist, NGenes, NBITS);

                    mutationFuncs[mutationChoice](children.second, generator, mutationDist, NGenes, NBITS);

                    size_t hash = children.first.computeHash ();
                    if(control.count(hash) == 0)
                    {
                        newPop.push_back(children.first);
                        control.insert(hash);
                    }
                    hash = children.second.computeHash();
                    if(control.count(hash) == 0)
                    {
                        newPop.push_back(children.second);
                        control.insert(hash);
                    }

                }
                if(newPop.size () >= popSize)
                {
                    break;
                }
                std::shuffle(selected.begin (), selected.end(), generator);
            }

            evaluate(newPop);
            pop.insert(pop.begin(), newPop.begin(), newPop.end());

            nondominatedSort();
            do
            {
                control.erase(pop.back().computeHash());
                pop.pop_back();
            } while(pop.size() > popSize);

            computeStatistics();
            computeCDF();
            onlinePerformance += stats.mean;
            offlinePerformance += stats.max;
        }
        onlinePerformance /= iterations;
        offlinePerformance /= iterations;

    }

    double getOnlinePerformance() const
    {
        return onlinePerformance;
    }

    double getOfflinePerformance() const
    {
        return offlinePerformance;
    }

    size_t getNumberOfSelectionOperators() const
    {
        return selectionFuncs.size();
    }

    size_t getNumberOfCrossoverOperators() const
    {
        return crossoverFuncs.size();
    }

    size_t getNumberOfMutationOperators() const
    {
        return mutationFuncs.size();
    }

    Statistics getStatistics() const
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
        return crossoverFuncNames[crossoverChoice];
    }
    string getMutationOperatorName()
    {
        return mutationFuncNames[mutationChoice];
    }

    vector<Chromosome> getBestChromosomes(size_t const n)
    {
        auto cmp = [&] (Individual const& lhs, Individual const& rhs)
        {
            if(lhs.rank < rhs.rank)
            {
                return true;
            }
            else if(lhs.rank == rhs.rank)
            {
                return lhs.distance > rhs.distance;
            }
            else
            {
                return false;
            }
        };

        size_t const sz = std::min(n, pop.size());
        std::partial_sort(pop.begin(), pop.begin() + sz, pop.end(), cmp);
        vector<Chromosome> res;
        res.reserve(sz);
        for(size_t i = 0; i < sz; ++i)
        {
            res.push_back(pop[i].chromosome);
        }
        return res;
    }

    vector<Individual> const& getPopulation() const
    {
        return pop;
    }

    vector<unsigned long> getDecodedChromosomes(size_t const n)
    {
        vector<unsigned long> res;
        size_t const sz = std::min(n, pop.size());
        res.reserve(sz*NGenes);
        for(size_t i = 0; i < sz; ++i)
        {
            for(size_t j = 0; j < NGenes; ++j)
            {
                res.push_back(pop[i].chromosome[j].to_ulong());
            }
        }
        return res;
    }

    size_t getNumberOfGenes() const
    {
        return NGenes;
    }




};


#endif // SOGA

