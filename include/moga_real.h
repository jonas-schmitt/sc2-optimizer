#ifndef _MOGA_REAL_
#define _MOGA_REAL_

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
#include <set>
#include <mpi.h>
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
using std::set;
using std::function;
using std::cout;
using std::endl;

struct SelectionParameter
{
    SelectionParameter(size_t N_, mt19937& generator_, vector<Individual>& pop_, size_t tournamentSize_)
        : N(N_), generator(&generator_), pop(&pop_),  tournamentSize(tournamentSize_){}
    size_t N;
    mt19937 *generator;
    vector<Individual> *pop;
    size_t tournamentSize;

};

struct CrossoverParameter
{
    CrossoverParameter(vector<Individual> const& parents_, mt19937& generator_, size_t NGenes_, size_t crossoverPoints_)
        : parents(&parents_), generator(&generator_), NGenes(NGenes_), crossoverPoints(crossoverPoints_) {}
    vector<Individual> const * parents;
    mt19937 *generator;
    size_t NGenes;
    size_t crossoverPoints;
};

struct MutationParameter
{
    MutationParameter(Individual& individual_, mt19937& generator_, std::bernoulli_distribution& mutationDist_, size_t NGenes_)
        : individual(&individual_), generator(&generator_), mutationDist(&mutationDist_), NGenes(NGenes_) {}
    Individual* individual;
    mt19937* generator;
    std::bernoulli_distribution * mutationDist;
    size_t NGenes;
};

template <typename T, typename U, Player const player>
class MOGA_real final
{
private:
    double mMutationProbability = 0.01;

    Statistics mStats;

    size_t mPopSize;

    size_t mNGenes;

    mt19937 mGenerator;
    bernoulli_distribution mFlipCoin;
    std::uniform_int_distribution<size_t> mChooseIndividual;
    std::uniform_real_distribution<double> mSpinWheel;
    bernoulli_distribution mMutationDist;

    vector<Individual> mPop;

    vector<vector<MicroSimulation<T,U>>> mSims;

    unordered_set<size_t> mPopControl;

    vector<string> mSelectionFuncNames = {"Tournament Selection", "Roulette Wheel Selection", "Stochastic Universal Sampling"};
    vector<string> mCrossoverFuncNames = {"N-Point Crossover", "Uniform Crossover"};
    vector<string> mMutationFuncNames = {"Uniform Mutation", "Gaussian Mutation"};

    size_t mSelectionChoice = 0;
    size_t mCrossoverChoice = 0;
    size_t mMutationChoice = 0;

    double mOnlinePerformance = 0.0;
    double mOfflinePerformance = 0.0;

    size_t mNGoals;

    size_t mTournamentSize = 4;
    size_t mNCrossoverPoints;



    function<vector<Individual *>(SelectionParameter)> tournamentSelection = [&](SelectionParameter params)
    {

        size_t const N = params.N;
        mt19937& generator = *params.generator;
        vector<Individual>& pop = *params.pop;
        size_t tournamentSize = params.tournamentSize;
        std::uniform_int_distribution<size_t> dist(0,pop.size()-1);
        auto func = [&] (std::uniform_int_distribution<size_t>& dist, mt19937& generator)
        {
            vector<size_t> positions;
            positions.reserve(tournamentSize);
            do
            {
                positions.push_back(dist(generator));
            } while(positions.size() < tournamentSize);
            while(positions.size() > 1)
            {
                Individual const& ind1 = pop[positions.back()];
                Individual const& ind2 = pop[positions[positions.size()-2]];
                if(ind1.rank < ind2.rank)
                {
                    positions[positions.size()-2] = positions.back();
                }
                else if(ind1.rank == ind2.rank && ind1.distance > ind2.distance)
                {
                    positions[positions.size()-2] = positions.back();
                }
                positions.pop_back();
            }
            return positions.front();

        };

        vector<Individual *> res;
        res.reserve(N);
        for(size_t i = 0; i < N; ++i)
        {
            res.push_back(&pop[func(dist, generator)]);
        }
        return res;
    };


    function<vector<Individual *>(SelectionParameter)> rouletteWheelSelection = [&](SelectionParameter params)
    {
        size_t const N = params.N;
        mt19937& generator = *params.generator;
        vector<Individual>& pop = *params.pop;
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

    function<vector<Individual *>(SelectionParameter)> stochasticUniversalSampling = [&](SelectionParameter params)
    {
        size_t const N = params.N;
        mt19937& generator = *params.generator;
        vector<Individual>& pop = *params.pop;
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



    function<pair<Individual, Individual>(CrossoverParameter)> nPointCrossover = [&] (CrossoverParameter params)
    {
        vector<Individual> const& parents = *params.parents;
        mt19937 &generator = *params.generator;
        size_t const NGenes = params.NGenes;
        size_t const crossoverPoints = params.crossoverPoints;

        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
        std::uniform_int_distribution<size_t> dist(1,NGenes-2);

        vector<size_t> posArray(crossoverPoints);
        unordered_set<size_t> positions;
        for(size_t& pos : posArray)
        {
            do {
            pos = dist(generator);
            } while(positions.count(pos) == 1);
            positions.insert(pos);
        }
        std::sort(posArray.begin(), posArray.end());
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
        for(size_t const pos : posArray)
        {
            for(size_t i = start; i < pos; ++i)
            {
                A.chromosome[i] = parent1.chromosome[i];
                B.chromosome[i] = parent2.chromosome[i];
            }
            swap(A,B);
            start = pos + 1;
        }
        for(size_t i = start; i < NGenes; ++i)
        {
            A.chromosome[i] = parent1.chromosome[i];
            B.chromosome[i] = parent2.chromosome[i];
        }
        return std::make_pair(child1, child2);
    };


    function<pair<Individual, Individual>(CrossoverParameter)> uniformCrossover = [&] (CrossoverParameter params)
    {
        vector<Individual> const& parents = *params.parents;
        mt19937 &generator = *params.generator;
        size_t const NGenes = params.NGenes;
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);

        std::bernoulli_distribution dist(0.5);
        Individual child1(NGenes), child2(NGenes);

        for(size_t i = 0; i < NGenes; ++i)
        {
            if(dist(generator))
            {
                child1.chromosome[i] = parent1.chromosome[i];
                child2.chromosome[i] = parent2.chromosome[i];
            }
            else
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
        }
        return std::make_pair(child1, child2);
    };



    function<void(MutationParameter)> uniformMutation = [&] (MutationParameter params)
    {
        Individual& individual = *params.individual;
        mt19937& generator = *params.generator;
        std::bernoulli_distribution& mutationDist = *params.mutationDist;
        size_t NGenes = params.NGenes;
        std::uniform_real_distribution<double> valueDist(MIN,MAX);
        for(size_t i = 0; i < NGenes; ++i)
        {
            if(mutationDist(generator))
            {
                individual.chromosome[i] = valueDist(generator);
            }
        }
    };

    function<void(MutationParameter)> gaussianMutation = [&] (MutationParameter params)
    {
        Individual& individual = *params.individual;
        mt19937& generator = *params.generator;
        std::bernoulli_distribution& mutationDist = *params.mutationDist;
        size_t NGenes = params.NGenes;
        std::normal_distribution<double> valueDist(0,0.5);
        for(size_t i = 0; i < NGenes; ++i)
        {
            if(mutationDist(generator))
            {
                double const val = individual.chromosome[i] + valueDist(generator);
                if(val < MIN) individual.chromosome[i] = MIN;
                else if(val > MAX) individual.chromosome[i] = MAX;
            }
        }
    };


    vector<function<vector<Individual *>(SelectionParameter) > > selectionFuncs = {tournamentSelection, rouletteWheelSelection, stochasticUniversalSampling};
    vector<function<pair<Individual, Individual>(CrossoverParameter)> > crossoverFuncs = {nPointCrossover, uniformCrossover};
    vector<function<void(MutationParameter) > > mutationFuncs = {uniformMutation, gaussianMutation};





    void evaluate(vector<Individual>& pop)
    {

        #pragma omp parallel for schedule(runtime)
        for(size_t i = 0; i < pop.size(); ++i)
        {
            double const value = 1.0/mNGoals;
            pop[i].fitness = 0;
            for(size_t j = 0; j < mNGoals; ++j)
            {
                if(player == Player::first)
                {
                    mSims[omp_get_thread_num()][j].setPlayer1Chromosome(pop[i].chromosome);
                    pop[i].fitness += mSims[omp_get_thread_num()][j].run(true, Player::first);
                }
                else
                {
                    mSims[omp_get_thread_num()][j].setPlayer2Chromosome(pop[i].chromosome);
                    pop[i].fitness += mSims[omp_get_thread_num()][j].run(true, Player::second);
                }

            }
            pop[i].fitness *= value;
        }
    }

    void setGoals(vector<Chromosome> const& goals)
    {
        if(goals.size() != mNGoals)
        {
            throw std::invalid_argument("MOGA_real::setGoals(): Invalid number of arguments");
        }
        #pragma omp parallel
        {
            for(size_t j = 0; j < mNGoals; ++j)
            {
                if(player == Player::first)
                {
                    mSims[omp_get_thread_num()][j].setPlayer2Chromosome(goals[j]);
                }
                else
                {
                    mSims[omp_get_thread_num()][j].setPlayer1Chromosome(goals[j]);
                }
            }
            #pragma omp barrier
        }
    }

    void computeCDF()
    {
        double const inv_sum = 1.0/mStats.sum;
        double tmp = 0.0;
        for(Individual& ind : mPop)
        {
            tmp += ind.fitness.score * inv_sum;
            ind.cdf = tmp;
            ind.total = mStats.sum;
        }
    }

    void computeStatistics()
    {
        vector<double> v(mPop.size());
        for(size_t i = 0; i < v.size(); ++i)
        {
            v[i] = mPop[i].fitness.score;
        }
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        double mean = sum / v.size();
        std::vector<double> diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(),
                       std::bind2nd(std::minus<double>(), mean));
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / v.size());
        mStats.sum = sum;
        mStats.mean = mean;
        mStats.max = *std::max_element(v.begin(), v.end());
        mStats.stdev = stdev;

    }

    void nondominatedSort()
    {
        // nondominated fronts
        vector<vector<size_t>> fronts(1);

        // determine the domination count and set of dominated individuals for all members of the population
        for(size_t i = 0; i < mPop.size(); ++i)
        {
            mPop[i].dominationCount = 0;
            mPop[i].dominationSet.clear();
            mPop[i].distance = 0.0;
            for(size_t j = 0; j < mPop.size(); ++j)
            {
                if(i == j) continue;
                if(mPop[i].dominates(mPop[j]))
                {
                    mPop[i].dominationSet.push_back(j);
                }
                else if(mPop[j].dominates(mPop[i]))
                {
                    ++mPop[i].dominationCount;
                }
            }
            // if the individual is not dominated by any other individuals, include it in the first nondominated front
            if(mPop[i].dominationCount == 0)
            {
                fronts[0].push_back(i);
                mPop[i].rank = 1;
            }
        }

        // determine the individuals in each subsequent front
        size_t r = 1;
        size_t count = fronts[0].size();
        fronts.emplace_back();
        for(size_t i = 0; count < mPopSize && !fronts[i].empty(); ++i)
        {
            for(size_t const p : fronts[i])
            {
                for(size_t const q : mPop[p].dominationSet)
                {
                    if(mPop[q].dominationCount > 0)
                    {
                        if(--mPop[q].dominationCount == 0)
                        {
                            fronts[i+1].push_back(q);
                            mPop[q].rank = r + 1;
                        }
                    }
                }
            }
            count += fronts[i+1].size();
            fronts.emplace_back();
            ++r;
        }


        vector<Individual> newPop;
        newPop.reserve(mPopSize);
        fronts.pop_back();


        auto cmp_damage = [&] (size_t const p, size_t const q)
        {
            return mPop[p].fitness.damage < mPop[q].fitness.damage;
        };
        auto cmp_health = [&] (size_t const p, size_t const q)
        {
            return mPop[p].fitness.health < mPop[q].fitness.health;
        };

        auto cmp_time = [&] (size_t const p, size_t const q)
        {
            return mPop[p].fitness.timeSteps < mPop[q].fitness.timeSteps;
        };


        // diversity preservation
        for(vector<size_t>& front : fronts)
        {

            // compute crowding distance for damage
            std::sort(front.begin(), front.end(), cmp_damage);
            mPop[front.front()].distance = INF;
            mPop[front.back()].distance = INF;
            if(front.size() > 2)
            {
                for(size_t i = 1; i < front.size()-2; ++i)
                {
                    mPop[front[i]].distance += mPop[front[i+1]].fitness.damage - mPop[front[i-1]].fitness.damage;
                }
            }



            // compute crowding distance for health
            std::sort(front.begin(), front.end(), cmp_health);
            mPop[front.front()].distance = INF;
            mPop[front.back()].distance = INF;
            if(front.size() > 2)
            {
                for(size_t i = 1; i < front.size()-2; ++i)
                {
                    if(mPop[front[i]].distance < INF)
                    {
                        mPop[front[i]].distance += mPop[front[i+1]].fitness.health - mPop[front[i-1]].fitness.health;
                    }
                }
            }

            // compute crowding distance for health
            std::sort(front.begin(), front.end(), cmp_time);
            mPop[front.front()].distance = INF;
            mPop[front.back()].distance = INF;
            if(front.size() > 2)
            {
                for(size_t i = 1; i < front.size()-2; ++i)
                {
                    if(mPop[front[i]].distance < INF)
                    {
                        mPop[front[i]].distance += mPop[front[i+1]].fitness.timeSteps - mPop[front[i-1]].fitness.timeSteps;
                    }
                }
            }


            // if the front does not fit completely in the population, sort the individuals according to the crowded comparison operator
            if(newPop.size() + front.size() > mPopSize)
            {
                auto cmp_distance = [&] (size_t const p, size_t const q)
                {
                    return mPop[p].distance > mPop[q].distance;
                };
                sort(front.begin(), front.end(), cmp_distance);
            }

            // consecutively include individuals from each front into the new population until it is filled
            for(size_t const p : front)
            {
                if(newPop.size() == mPopSize) break;
                newPop.push_back(mPop[p]);
            }

        }
        // replace the old population with the new one
        mPop = newPop;
    }

    size_t countUnitTypes(vector<string> buildOrder)
    {
        set<string> strSet(buildOrder.begin(), buildOrder.end());
        return strSet.size();
    }


public:

    typedef T race1;
    typedef U race2;

    MOGA_real(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2, size_t popSize, vector<string> const & buildList1, vector<string> const & buildList2, size_t const nGoals)
        :  mPopSize(popSize), mGenerator(std::chrono::system_clock::now().time_since_epoch().count()), mFlipCoin(0.5), mChooseIndividual(0,popSize-1), mSpinWheel(0,1.0), mNGoals(nGoals)
    {

        mSims.reserve(omp_get_num_threads());
        #pragma omp parallel
        {
            for(size_t i = 0; i < omp_get_num_threads(); ++i)
            {
                #pragma omp critical
                {
                    if(i == omp_get_thread_num())
                    {
                        mSims.emplace_back();
                        mSims[omp_get_thread_num()].reserve(mNGoals);
                        for(size_t j = 0; j < mNGoals; ++j)
                        {
                            mSims[omp_get_thread_num()].emplace_back(minPos, maxPos, filePath1, filePath2);
                            mSims[omp_get_thread_num()].back().initBothPlayers(buildList1, buildList2);
                        }
                    }
                }
                #pragma omp barrier
            }
        }

        if(player == Player::first)
        {
            mNGenes = mSims[0][0].getPlayer1ChromosomeLength();
            mNCrossoverPoints = countUnitTypes(buildList1);
        }
        else
        {
            mNGenes = mSims[0][0].getPlayer2ChromosomeLength();
            mNCrossoverPoints = countUnitTypes(buildList2);
        }
        if(mNCrossoverPoints > 1) --mNCrossoverPoints;
        std::uniform_real_distribution<double> valueDist(MIN, MAX);
        vector<Chromosome> initChroms(mNGoals);
        for(auto& chrom : initChroms)
        {
            chrom.resize(mNGenes);
            for(auto& gene : chrom)
            {
                gene = valueDist(mGenerator);
            }
        }
        setGoals(initChroms);


        mPop.reserve(2*popSize);
        mPop.resize(popSize);


        mStats.sum = 0.0;
        mStats.max = 0.0;

        for(size_t i = 0; i < popSize; ++i)
        {
            Chromosome& chrom = mPop[i].chromosome;
            chrom.resize(mNGenes);
            double hash;
            do
            {
                for(auto& gene : chrom)
                {
                    gene = valueDist(mGenerator);
                }
                hash = mPop[i].computeHash();
            }
            while(mPopControl.count(hash) == 1);
            mPopControl.insert(hash);
        }
        evaluate(mPop);
        computeStatistics();
        computeCDF();

    }

    void optimize(vector<Chromosome> const& goals, size_t const iterations)
    {
        mOnlinePerformance = 0.0;
        mOfflinePerformance = 0.0;

        mMutationDist = std::move(std::bernoulli_distribution(mMutationProbability));
        setGoals(goals);
        evaluate(mPop);

        computeStatistics();
        computeCDF();
        for(size_t i = 0; i < iterations; ++i)
        {
            // apply genetic algorithm

            vector<Individual> newPop;
            newPop.reserve(mPopSize);
            vector<Individual *> selected;

            selected = selectionFuncs[mSelectionChoice](SelectionParameter(mPopSize, mGenerator, mPop, mTournamentSize));

            for(size_t count = 0; count < 100; ++count)
            {
                if(selected.size() == 0) break; // Just for safety
                for(size_t i = 0; i < selected.size()-1 && newPop.size() < mPopSize; i += 2)
                {
                    pair<Individual, Individual> children;

                    children = crossoverFuncs[mCrossoverChoice](CrossoverParameter({*selected[i], *selected[i+1]}, mGenerator, mNGenes, mNCrossoverPoints));
                    mutationFuncs[mMutationChoice](MutationParameter(children.first, mGenerator, mMutationDist, mNGenes));

                    mutationFuncs[mMutationChoice](MutationParameter(children.second, mGenerator, mMutationDist, mNGenes));

                    size_t hash = children.first.computeHash ();
                    if(mPopControl.count(hash) == 0)
                    {
                        newPop.push_back(children.first);
                        mPopControl.insert(hash);
                    }
                    hash = children.second.computeHash();
                    if(mPopControl.count(hash) == 0)
                    {
                        newPop.push_back(children.second);
                        mPopControl.insert(hash);
                    }

                }
                if(newPop.size () >= mPopSize)
                {
                    break;
                }
                std::shuffle(selected.begin (), selected.end(), mGenerator);
            }

            evaluate(newPop);
            mPop.insert(mPop.begin(), newPop.begin(), newPop.end());

            nondominatedSort();
            do
            {
                mPopControl.erase(mPop.back().computeHash());
                mPop.pop_back();
            } while(mPop.size() > mPopSize);

            computeStatistics();
            computeCDF();
            mOnlinePerformance += mStats.mean;
            mOfflinePerformance += mStats.max;
        }
        mOnlinePerformance /= iterations;
        mOfflinePerformance /= iterations;

    }

    double getOnlinePerformance() const
    {
        return mOnlinePerformance;
    }

    double getOfflinePerformance() const
    {
        return mOfflinePerformance;
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
        return mStats;
    }

    void setSelection(size_t const value)
    {
        mSelectionChoice = value < selectionFuncs.size() ? value : mSelectionChoice;
    }

    void setCrossover(size_t const value)
    {
        mCrossoverChoice = value < crossoverFuncs.size() ? value : mCrossoverChoice;
    }

    void setMutation(size_t const value)
    {
        mMutationChoice = value < mutationFuncs.size() ? value : mMutationChoice;
    }

    string getSelectionOperatorName()
    {
        return mSelectionFuncNames[mSelectionChoice];
    }
    string getCrossoverOperatorName()
    {
        return mCrossoverFuncNames[mCrossoverChoice];
    }
    string getMutationOperatorName()
    {
        return mMutationFuncNames[mMutationChoice];
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

        size_t const sz = std::min(n, mPop.size());
        std::partial_sort(mPop.begin(), mPop.begin() + sz, mPop.end(), cmp);
        vector<Chromosome> res;
        res.reserve(sz);
        for(size_t i = 0; i < sz; ++i)
        {
            res.push_back(mPop[i].chromosome);
        }
        return res;
    }

    vector<Individual> const& getPopulation() const
    {
        return mPop;
    }

    size_t getPopulationSize() const
    {
        return mPop.size();
    }

    Chromosome getChromosomes(size_t const migrants)
    {
        Chromosome res;
        res.reserve(migrants*mNGenes);
        for(size_t i = 0; i < migrants; ++i)
        {
            for(size_t j = 0; j < mNGenes; ++j)
            {
                res.push_back(mPop[i].chromosome[j]);
            }
        }
        return res;
    }

    size_t getNumberOfGenes() const
    {
        return mNGenes;
    }

    void includeDecodedChromosomes(Chromosome const& data, size_t const migrants, int const rank, int const procs)
    {
        vector<Individual> newPop;
        newPop.reserve((procs-1)*migrants);
        for(int i = 0; i < procs; ++i)
        {
            if(rank == i) continue;
            for(size_t j = 0; j < migrants; ++j)
            {
                newPop.emplace_back();
                for(size_t k = 0; k < mNGenes; ++k)
                {
                    newPop.back().chromosome.emplace_back(data[i*migrants*mNGenes + j*mNGenes + k]);
                }

            }
        }
        evaluate(newPop);
        mPop.insert(mPop.begin(), newPop.begin(), newPop.end());

        nondominatedSort();
        do
        {
            mPopControl.erase(mPop.back().computeHash());
            mPop.pop_back();
        } while(mPop.size() > mPopSize);
    }




};


#endif // SOGA

