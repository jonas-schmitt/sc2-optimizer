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
#include <set>
#include <mpi.h>
#include <omp.h>
#include <functional>


#include "Chromosome.h"
#include "MicroSimulation.h"
#include "Utilities.h"


struct CrossoverParameter
{
    CrossoverParameter(std::vector<Individual> const& parents_, std::mt19937& generator_, size_t crossoverPoints_)
        : parents(&parents_), generator(&generator_), crossoverPoints(crossoverPoints_) {}
    std::vector<Individual> const * parents;
    std::mt19937 *generator;
    size_t crossoverPoints;
};

struct MutationParameter
{
    MutationParameter(Individual& individual_, std::mt19937& generator_, size_t geneToMutate_, size_t currentGeneration_, size_t maxGenerations_)
        : individual(&individual_), generator(&generator_), geneToMutate(geneToMutate_), currentGeneration(currentGeneration_), maxGenerations(maxGenerations_){}
    Individual* individual;
    std::mt19937* generator;
    size_t geneToMutate;
    size_t currentGeneration;
    size_t maxGenerations;
};

template <typename T, typename U, Player const player>
class MOGA final
{
private:
    double mMutationProbability;

    std::pair<Statistics,Statistics> mStats;

    size_t mPopSize;

    size_t mNGenes;


    std::uniform_real_distribution<double> mDistribution;

    std::vector<Individual> mPop;

    std::vector<std::vector<MicroSimulation<T,U>>> mSims;



    //std::vector<std::string> mSelectionFuncNames = {"Tournament Selection", "Roulette Wheel Selection", "Stochastic Universal Sampling"};
    std::vector<std::string> mCrossoverFuncNames = {"Self-Adaptive Simulated Binary Crossover", "Simulated Binary Crossover", "N-Point Crossover", "Uniform Crossover", "Intermediate Crossover", "Line Crossover", "Arithmetic Crossover"};
    std::vector<std::string> mMutationFuncNames = {"Polynomial Mutation", "Gaussian Mutation", "Uniform Mutation"};

    //size_t mSelectionChoice = 0;
    size_t mCrossoverChoice = 0;
    size_t mMutationChoice = 0;

    size_t mNGoals;

    size_t mTournamentSize = 4;
    size_t mNCrossoverPoints;

    size_t mIndividualToMutate;
    size_t mGeneToMutate;

    size_t mMutationCount = 0;

    // Crossover methods


    std::function<std::pair<Individual, Individual>(CrossoverParameter)> simulatedBinaryCrossover = [&] (CrossoverParameter params)
    {
        std::vector<Individual> const& parents = *params.parents;
        std::mt19937 &generator = *params.generator;
        std::uniform_real_distribution<double>& dist = mDistribution;
        size_t const NGenes = mNGenes;
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);

        Individual child1(NGenes), child2(NGenes);


        double constexpr n_c = 2.0;
        double constexpr exp1 = n_c + 1.0;
        double constexpr exp2 = 1.0 / exp1;
        double const u = dist(generator);
        double const a = u < 0.5 ? std::pow(2.0 * u, exp2) : std::pow(1.0 / (2.0 * (1.0 - u)), exp2);

        for(size_t i = 0; i < NGenes; ++i)
        {
            double const x1 = parent1.chromosome[i];
            double const x2 = parent2.chromosome[i];

            double const avg = 0.5 * (x1 + x2);
            double const diff = 0.5 * (x2 - x1);
            double const y1 = avg - a * diff;
            double const y2 = avg + a * diff;
            child1.chromosome[i] = y1;
            child2.chromosome[i] = y2;

        }
        return std::make_pair(child1, child2);

    };

    std::function<std::pair<Individual, Individual>(CrossoverParameter)> adaptiveSBX = [&] (CrossoverParameter params)
    {
        std::vector<Individual> const& parents = *params.parents;
        std::mt19937 &generator = *params.generator;
        std::uniform_real_distribution<double>& dist = mDistribution;
        size_t const NGenes = mNGenes;
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);

        Individual child1(NGenes), child2(NGenes);
        double const u = dist(generator);

        bool const contraction = u < 0.5;

        double constexpr n_c = 2.0;
        double constexpr exp = 1.0/(n_c + 1.0);

        double const tmp1 = 2.0 * u;
        double const tmp2 = 1.0 / (2.0 * (1.0 - u));

        double const a = contraction ? std::pow(tmp1, exp) : std::pow(tmp2, exp);

        for(size_t i = 0; i < NGenes; ++i)
        {
            double const x1 = parent1.chromosome[i];
            double const x2 = parent2.chromosome[i];

            double const avg = 0.5 * (x1 + x2);
            double const diff = 0.5 * (x2 - x1);
            double const y1 = avg - a * diff;
            double const y2 = avg + a * diff;
            child1.chromosome[i] = y1;
            child2.chromosome[i] = y2;

        }

        child1.fitness = 0;
        child2.fitness = 0;
        for(size_t j = 0; j < mNGoals; ++j)
        {
            if(player == Player::first)
            {
                mSims[omp_get_thread_num()][j].setPlayer1Chromosome(child1.chromosome);
                child1.fitness += mSims[omp_get_thread_num()][j].run(true, Player::first);
                mSims[omp_get_thread_num()][j].setPlayer1Chromosome(child2.chromosome);
                child2.fitness += mSims[omp_get_thread_num()][j].run(true, Player::first);
            }
            else
            {
                mSims[omp_get_thread_num()][j].setPlayer2Chromosome(child1.chromosome);
                child1.fitness += mSims[omp_get_thread_num()][j].run(true, Player::second);
                mSims[omp_get_thread_num()][j].setPlayer2Chromosome(child2.chromosome);
                child2.fitness += mSims[omp_get_thread_num()][j].run(true, Player::second);
            }
        }
        child1.fitness /= static_cast<double>(mNGoals);
        child2.fitness /= static_cast<double>(mNGoals);
        child1.evaluated = true;
        child1.evaluated = true;

        // First child


        bool const worse1 = parent1.dominates(child1) && parent2.dominates(child1);
        bool const better1 = child1.dominates(parent1) && child1.dominates(parent2);
        bool const worse2 = parent1.dominates(child2) && parent2.dominates(child2);
        bool const better2 = child2.dominates(parent1) && child2.dominates(parent2);

        if(better1 || worse1)
        {
            child1.evaluated = false;
        }
        if(better2 || worse2)
        {
            child2.evaluated = false;
        }

        double constexpr alpha = 1.5;

        bool skip1 = false, skip2 = false;
        if(contraction)
        {
            double exp1, exp2, a1, a2;
            if(better1)
            {
                double const n_c1 = std::max(0.0, std::min(50.0, (1.0+n_c)/alpha - 1.0));
                exp1 = 1.0 / (n_c1 + 1.0);
                a1 = std::pow(tmp1, exp1);
            }
            else if(worse1)
            {
                double const n_c1 = std::max(0.0, std::min(50.0, alpha*(1.0+n_c) - 1.0));
                exp1 = 1.0 / (n_c1 + 1.0);
                a1 = std::pow(tmp1, exp1);
            }
            else
            {
                skip1 = true;
            }

            if(better2)
            {
                if(better1)
                {
                    exp2 = exp1;
                    a2 = a1;
                }
                else
                {
                    double const n_c2 = std::max(0.0, std::min(50.0, (1.0+n_c)/alpha - 1.0));
                    exp2 = 1.0 / (n_c2 + 1.0);
                    a2 = std::pow(tmp1, exp2);
                }
            }
            else if(worse2)
            {
                if(worse1)
                {
                    exp2 = exp1;
                    a2 = a1;
                }

                else
                {
                    double const n_c2 = std::max(0.0, std::min(50.0, alpha*(1.0+n_c) - 1.0));
                    exp2 = 1.0 / (n_c2 + 1.0);
                    a2 = std::pow(tmp1, exp2);
                }
            }
            else
            {
               skip2 = true;
            }

            if(!skip1 && !skip2)
            {
                child1.alternative = child1.chromosome;
                child2.alternative = child2.chromosome;
                for(size_t i = 0; i < NGenes; ++i)
                {
                    double const x1 = parent1.chromosome[i];
                    double const x2 = parent2.chromosome[i];

                    double const avg = 0.5 * (x1 + x2);
                    double const diff = 0.5 * (x2 - x1);
                    double const y1 = avg - a1 * diff;
                    double const y2 = avg + a2 * diff;
                    child1.chromosome[i] = y1;
                    child2.chromosome[i] = y2;
                }
            }
            else if(skip1)
            {
                child2.alternative = child2.chromosome;
                for(size_t i = 0; i < NGenes; ++i)
                {
                    double const x1 = parent1.chromosome[i];
                    double const x2 = parent2.chromosome[i];

                    double const avg = 0.5 * (x1 + x2);
                    double const diff = 0.5 * (x2 - x1);
                    double const y2 = avg + a2 * diff;
                    child2.chromosome[i] = y2;

                }
            }
            else if(skip2)
            {
                child1.alternative = child1.chromosome;
                for(size_t i = 0; i < NGenes; ++i)
                {
                    double const x1 = parent1.chromosome[i];
                    double const x2 = parent2.chromosome[i];

                    double const avg = 0.5 * (x1 + x2);
                    double const diff = 0.5 * (x2 - x1);
                    double const y1 = avg - a1 * diff;
                    child1.chromosome[i] = y1;
                }
            }
        }
        else
        {

            if(better1)
            {
                if(better2)
                {
                    child1.alternative = child1.chromosome;
                    child2.alternative = child2.chromosome;
                    for(size_t i = 0; i < NGenes; ++i)
                    {
                        double const x1 = parent1.chromosome[i];
                        double const x2 = parent2.chromosome[i];
                        double const diff = x2 - x1;
                        double const beta = std::abs((child2.chromosome[i] - child1.chromosome[i])/(diff));
                        double n_c12 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * std::log(beta)/std::log(1 + alpha*(beta - 1))));


                        double const exp = 1.0 / (n_c12 + 1.0);
                        double const a = std::pow(tmp2, exp);


                        double const avg = 0.5 * (x1 + x2);
                        double const shift = a * 0.5 * diff;
                        double const y1 = avg - shift;
                        double const y2 = avg + shift;
                        child1.chromosome[i] = y1;
                        child2.chromosome[i] = y2;
                    }
                }
                else if(worse2)
                {
                    child1.alternative = child1.chromosome;
                    child2.alternative = child2.chromosome;
                    for(size_t i = 0; i < NGenes; ++i)
                    {
                        double const x1 = parent1.chromosome[i];
                        double const x2 = parent2.chromosome[i];
                        double const diff = x2 - x1;
                        double const beta = std::abs((child2.chromosome[i] - child1.chromosome[i])/(diff));

                        double const log_beta = std::log(beta);
                        double constexpr alpha_inv = 1.0/alpha;
                        double const n_c1 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * log_beta/std::log(1 + alpha*(beta - 1))));

                        double const n_c2 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * log_beta/std::log(1 + alpha_inv*(beta - 1))));

                        double const exp1 = 1.0 / (n_c1 + 1.0);
                        double const a1 = std::pow(tmp2, exp1);

                        double const exp2 = 1.0 / (n_c2 + 1.0);
                        double const a2 = std::pow(tmp2, exp2);

                        double const avg = 0.5 * (x1 + x2);
                        double const y1 = avg - a1 * 0.5 * diff;
                        double const y2 = avg + a2 * 0.5 * diff;
                        child1.chromosome[i] = y1;
                        child2.chromosome[i] = y2;

                    }
                }
                else
                {
                    child1.alternative = child1.chromosome;
                    for(size_t i = 0; i < NGenes; ++i)
                    {
                        double const x1 = parent1.chromosome[i];
                        double const x2 = parent2.chromosome[i];
                        double const diff = x2 - x1;
                        double const beta = std::abs((child2.chromosome[i] - child1.chromosome[i])/(diff));
                        double n_c1 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * std::log(beta)/std::log(1 + alpha*(beta - 1))));

                        double const exp = 1.0 / (n_c1 + 1.0);
                        double const a = std::pow(tmp2, exp);

                        double const avg = 0.5 * (x1 + x2);

                        double const y1 = avg - a * 0.5 * diff;
                        child1.chromosome[i] = y1;
                    }
                }
            }
            else if(worse1)
            {
                if(better2)
                {
                    child1.alternative = child1.chromosome;
                    child2.alternative = child2.chromosome;
                    for(size_t i = 0; i < NGenes; ++i)
                    {
                        double const x1 = parent1.chromosome[i];
                        double const x2 = parent2.chromosome[i];
                        double const diff = x2 - x1;
                        double const beta = std::abs((child2.chromosome[i] - child1.chromosome[i])/(diff));

                        double const log_beta = std::log(beta);
                        double constexpr alpha_inv = 1.0/alpha;

                        double const n_c1 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * log_beta/std::log(1 + alpha_inv*(beta - 1))));

                        double const n_c2 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * log_beta/std::log(1 + alpha*(beta - 1))));

                        double const exp1 = 1.0 / (n_c1 + 1.0);
                        double const a1 = std::pow(tmp2, exp1);

                        double const exp2 = 1.0 / (n_c2 + 1.0);
                        double const a2 = std::pow(tmp2, exp2);


                        double const avg = 0.5 * (x1 + x2);
                        double const y1 = avg - a1 * 0.5 * diff;
                        double const y2 = avg + a2 * 0.5 * diff;
                        child1.chromosome[i] = y1;
                        child2.chromosome[i] = y2;

                    }
                }
                else if(worse2)
                {
                    child1.alternative = child1.chromosome;
                    child2.alternative = child2.chromosome;
                    for(size_t i = 0; i < NGenes; ++i)
                    {
                        double const x1 = parent1.chromosome[i];
                        double const x2 = parent2.chromosome[i];
                        double const diff = x2 - x1;
                        double const beta = std::abs((child2.chromosome[i] - child1.chromosome[i])/(diff));
                        double const n_c12 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * std::log(beta)/std::log(1 + (beta - 1)/alpha)));

                        double const exp = 1.0 / (n_c12 + 1.0);
                        double const a = std::pow(tmp2, exp);


                        double const avg = 0.5 * (x1 + x2);
                        double const shift = a * 0.5 * diff;
                        double const y1 = avg - shift;
                        double const y2 = avg + shift;
                        child1.chromosome[i] = y1;
                        child2.chromosome[i] = y2;

                    }
                }
                else
                {
                    child1.alternative = child1.chromosome;
                    for(size_t i = 0; i < NGenes; ++i)
                    {
                        double const x1 = parent1.chromosome[i];
                        double const x2 = parent2.chromosome[i];
                        double const diff = x2 - x1;
                        double const beta = std::abs((child2.chromosome[i] - child1.chromosome[i])/(diff));
                        double const n_c1 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * std::log(beta)/std::log(1 + (beta - 1)/alpha)));

                        double const exp = 1.0 / (n_c1 + 1.0);
                        double const a = std::pow(tmp2, exp);

                        double const avg = 0.5 * (x1 + x2);
                        double const y1 = avg - a * 0.5 * diff;
                        child1.chromosome[i] = y1;
                    }
                }
            }
            else
            {
                if(better2)
                {
                    child2.alternative = child2.chromosome;
                    for(size_t i = 0; i < NGenes; ++i)
                    {
                        double const x1 = parent1.chromosome[i];
                        double const x2 = parent2.chromosome[i];
                        double const diff = x2 - x1;
                        double const beta = std::abs((child2.chromosome[i] - child1.chromosome[i])/(diff));

                        double const n_c2 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * std::log(beta)/std::log(1 + alpha*(beta - 1))));

                        double const exp = 1.0 / (n_c2 + 1.0);
                        double const a = std::pow(tmp2, exp);

                        double const avg = 0.5 * (x1 + x2);
                        double const y2 = avg + a * 0.5* diff;
                        child2.chromosome[i] = y2;
                    }
                }
                else if(worse2)
                {
                    child2.alternative = child2.chromosome;
                    for(size_t i = 0; i < NGenes; ++i)
                    {
                        double const x1 = parent1.chromosome[i];
                        double const x2 = parent2.chromosome[i];
                        double const diff = x2 - x1;
                        double const beta = std::abs((child2.chromosome[i] - child1.chromosome[i])/(diff));

                        double const n_c2 = std::max(0.0, std::min(50.0, -1 + (n_c + 1) * std::log(beta)/std::log(1 + (beta - 1)/alpha)));

                        double const exp = 1.0 / (n_c2 + 1.0);
                        double const a = std::pow(tmp2, exp);

                        double const avg = 0.5 * (x1 + x2);
                        double const y2 = avg + 0.5 * a * diff;
                        child2.chromosome[i] = y2;
                    }
                }
            }
        }

        return std::make_pair(child1, child2);


    };


    std::function<std::pair<Individual, Individual>(CrossoverParameter)> nPointCrossover = [&] (CrossoverParameter params)
    {
        std::vector<Individual> const& parents = *params.parents;
        std::mt19937 &generator = *params.generator;
        size_t const NGenes = mNGenes;
        size_t const crossoverPoints = params.crossoverPoints;

        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);
        std::uniform_int_distribution<size_t> dist(1,NGenes-2);

        std::vector<size_t> posArray(crossoverPoints);
        std::unordered_set<size_t> positions;
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


    std::function<std::pair<Individual, Individual>(CrossoverParameter)> uniformCrossover = [&] (CrossoverParameter params)
    {
        std::vector<Individual> const& parents = *params.parents;
        std::mt19937 &generator = *params.generator;
        size_t const NGenes = mNGenes;
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

    std::function<std::pair<Individual, Individual>(CrossoverParameter)> intermediateCrossover = [&] (CrossoverParameter params)
    {
        std::vector<Individual> const& parents = *params.parents;
        std::mt19937 &generator = *params.generator;
        size_t const NGenes = mNGenes;
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);

        std::uniform_real_distribution<double> dist(0.25,1.25);
        Individual child1(NGenes), child2(NGenes);

        for(size_t i = 0; i < NGenes; ++i)
        {
            child1.chromosome[i] = parent1.chromosome[i] + dist(generator)*(parent2.chromosome[i] - parent1.chromosome[i]);
            child2.chromosome[i] = parent1.chromosome[i] + dist(generator)*(parent2.chromosome[i] - parent1.chromosome[i]);
        }
        return std::make_pair(child1, child2);
    };

    std::function<std::pair<Individual, Individual>(CrossoverParameter)> lineCrossover = [&] (CrossoverParameter params)
    {
        std::vector<Individual> const& parents = *params.parents;
        std::mt19937 &generator = *params.generator;
        size_t const NGenes = mNGenes;
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);

        std::uniform_real_distribution<double> dist(0.25,1.25);
        Individual child1(NGenes), child2(NGenes);

        double const alpha1 = dist(generator);
        double const alpha2 = dist(generator);

        for(size_t i = 0; i < NGenes; ++i)
        {
            child1.chromosome[i] = parent1.chromosome[i] + alpha1*(parent2.chromosome[i] - parent1.chromosome[i]);
            child2.chromosome[i] = parent1.chromosome[i] + alpha2*(parent2.chromosome[i] - parent1.chromosome[i]);
        }
        return std::make_pair(child1, child2);
    };

    std::function<std::pair<Individual, Individual>(CrossoverParameter)> arithmeticCrossover = [&] (CrossoverParameter params)
    {
        std::vector<Individual> const& parents = *params.parents;
        std::mt19937 &generator = *params.generator;
        std::uniform_real_distribution<double>& dist = mDistribution;
        size_t const NGenes = mNGenes;
        Individual const& parent1 = parents.at(0);
        Individual const& parent2 = parents.at(1);

        Individual child1(NGenes), child2(NGenes);

        double const alpha1 = dist(generator);
        double const alpha2 = dist(generator);

        for(size_t i = 0; i < NGenes; ++i)
        {
            child1.chromosome[i] = alpha1*parent1.chromosome[i] + (1.0-alpha1)*parent2.chromosome[i];
            child2.chromosome[i] = alpha2*parent1.chromosome[i] + (1.0-alpha2)*parent2.chromosome[i];
        }
        return std::make_pair(child1, child2);
    };


    std::function<void(MutationParameter)> polynomialMutation = [&] (MutationParameter params)
    {
        Individual& individual = *params.individual;
        std::mt19937& generator = *params.generator;
        std::uniform_real_distribution<double>& dist = mDistribution;
        size_t const i = params.geneToMutate;

        // comment in for adaptive mutation probability
        //double const t = static_cast<double>(params.currentGeneration);
        //double const t_max = static_cast<double>(params.maxGenerations);

        //double const len_inv = 1.0 / static_cast<double>(individual.chromosome.size());
        //double const prob = len_inv + (t / t_max)*(1 - len_inv);
        //std::bernoulli_distribution mutationDist(prob);

        double constexpr n_m = 100.0; // + t;
        double constexpr exp1 = n_m + 1;
        double constexpr exp2 = 1/exp1;
        double constexpr max_range = MAX - MIN;
        double constexpr max_range_inv = 1.0/max_range;

        double const u = dist(generator);

        double const a = std::min(individual.chromosome[i] - MIN, MAX - individual.chromosome[i]) * max_range_inv;
        double const b = std::pow(1.0 - a, exp1);

        double delta;
        if(u < 0.5)
        {
            double const c = 2.0 * u;
            double const d = 1.0 - c;
            delta = std::pow(c + d * b, exp2) - 1.0;
        }
        else
        {
            double const c = 2.0*(1.0 - u);
            double const d = 2.0*(u - 0.5);
            delta = 1.0 - std::pow(c + d * b, exp2);
        }
        individual.chromosome[i] += delta * max_range;

    };

    std::function<void(MutationParameter)> gaussianMutation = [&] (MutationParameter params)
    {
        Individual& individual = *params.individual;
        std::mt19937& generator = *params.generator;
        std::uniform_real_distribution<double>& dist = mDistribution;
        size_t const i = params.geneToMutate;
        double constexpr sigma = 0.375/(MAX - MIN);
        double constexpr tmp = std::sqrt(2.0)*(MAX-MIN)*sigma;


        double const u = dist(generator);
        double const x = individual.chromosome[i];
        double const u_min = 0.5 * std::erf((MIN - x)/tmp) + 0.5;
        double const u_max = 0.5 * std::erf((MAX - x)/tmp) + 0.5;
        double const a = u < 0.5 ? 2*u_min*(1 - 2*u) : 2*u_max*(2*u - 1);
        double constexpr b = 8*(M_PI - 3)/(3*M_PI*(4 - M_PI));
        double const c = std::log(1-a*a);
        double constexpr d = 2/(b*M_PI);
        double const e = d + c/2;
        double const f = c/b;
        double const g = sgn(a) * std::sqrt(std::sqrt(e*e - f) - e);

        individual.chromosome[i] += std::sqrt(2) * sigma * (MAX - MIN) * g;
    };


    std::function<void(MutationParameter)> uniformMutation = [&] (MutationParameter params)
    {
        Individual& individual = *params.individual;
        std::mt19937& generator = *params.generator;
        size_t const i = params.geneToMutate;
        std::uniform_real_distribution<double> valueDist(MIN,MAX);
        individual.chromosome[i] = valueDist(generator);
    };


    std::vector<std::function<std::pair<Individual, Individual>(CrossoverParameter)> > crossoverFuncs = {adaptiveSBX, simulatedBinaryCrossover, nPointCrossover, uniformCrossover, intermediateCrossover, lineCrossover, arithmeticCrossover};
    std::vector<std::function<void(MutationParameter) > > mutationFuncs = {polynomialMutation, gaussianMutation, uniformMutation};

    void mutationClock(std::mt19937& generator)
    {
        double const u = mDistribution(generator);
        double const l = mNGenes * (-std::log(1-u));
        size_t const tmp = static_cast<size_t>(mGeneToMutate + l);
        mIndividualToMutate = tmp / mNGenes;
        mGeneToMutate = tmp % mNGenes;
    }


    std::vector<Individual *> binaryTournamentSelection(size_t const N, std::mt19937& generator)
    {
        std::uniform_int_distribution<size_t> dist(0, mPop.size()-1);
        std::vector<Individual *> res;
        res.reserve(N);

        #pragma omp parallel
        {
            std::vector<Individual *> res_local;
            res_local.reserve(N / omp_get_num_threads() + N % omp_get_max_threads());

            #pragma omp for schedule(static)
            for(size_t i = 0; i < N; ++i)
            {

                size_t const pos1 = dist(generator);
                size_t pos2;
                while((pos2 = dist(generator)) == pos1);
                Individual& ind1 = mPop[pos1];
                Individual& ind2 = mPop[pos2];
                if(ind1.rank < ind2.rank || ind1.rank == ind2.rank && ind1.distance > ind2.distance)
                {
                    res_local.push_back(&ind1);
                }
                else
                {
                    res_local.push_back(&ind2);
                }
            }
            #pragma omp critical
            {
                res.insert(res.end(), res_local.begin(), res_local.end());
            }
        }
        return res;
    }

    std::vector<Individual *> tournamentSelection(size_t const N, std::mt19937& generator)
    {
        std::uniform_int_distribution<size_t> dist(0,mPop.size()-1);
        std::vector<Individual *> res;
        res.reserve(N);
        #pragma omp parallel
        {
            std::vector<Individual *> res_local;
            res_local.reserve(N / omp_get_num_threads() + N % omp_get_max_threads());
            std::vector<size_t> positions;
            positions.reserve(mTournamentSize);
            #pragma omp for schedule(static)
            for(size_t i = 0; i < N; ++i)
            {
                do
                {
                    positions.push_back(dist(generator));
                } while(positions.size() < mTournamentSize);
                while(positions.size() > 1)
                {
                    Individual const& ind1 = mPop[positions.back()];
                    Individual const& ind2 = mPop[positions[positions.size()-2]];
                    if(ind1.rank < ind2.rank || ind1.rank == ind2.rank && ind1.distance > ind2.distance)
                    {
                        positions[positions.size()-2] = positions.back();
                    }
                    positions.pop_back();
                }
                res_local.push_back(&mPop[positions[0]]);
                positions.pop_back();
            }
            #pragma omp critical
            {
                res.insert(res.end(), res_local.begin(), res_local.end());
            }
        }
        return res;
    }

    void evaluate(std::vector<Individual>& pop)
    {


        #pragma omp parallel for schedule(dynamic,1)
        for(size_t i = 0; i < pop.size(); ++i)
        {
            if(!pop[i].evaluated)
            {
                Fitness fn;
                for(size_t j = 0; j < mNGoals; ++j)
                {
                    if(player == Player::first)
                    {
                        mSims[omp_get_thread_num()][j].setPlayer1Chromosome(pop[i].chromosome);
                        fn += mSims[omp_get_thread_num()][j].run(true, Player::first);
                    }
                    else
                    {
                        mSims[omp_get_thread_num()][j].setPlayer2Chromosome(pop[i].chromosome);
                        fn += mSims[omp_get_thread_num()][j].run(true, Player::second);
                    }

                }
                fn /= static_cast<double>(mNGoals);
                if(!pop[i].alternative.empty())
                {
                    if(pop[i].fitness.dominates(fn))
                    {
                        pop[i].chromosome = pop[i].alternative;
                    }
                    else
                    {
                        pop[i].fitness = fn;
                    }
                    pop[i].alternative.clear();
                }
		else
		{
		    pop[i].fitness = fn;
		}
            }
        }
    }

    void setGoals(std::vector<Chromosome> const& goals)
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
        }
    }


    void computeStatistics(Statistics& stats, std::vector<double> const& v)
    {
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
        stats.onlinePerformance += stats.mean;
        stats.offlinePerformance += stats.max;
        ++stats.iteration;

    }

    void nondominatedSort()
    {
        // nondominated fronts
        std::vector<std::vector<size_t>> fronts(1);

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


        std::vector<Individual> newPop;
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
        for(std::vector<size_t>& front : fronts)
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

    size_t countUnitTypes(std::vector<std::string> buildOrder)
    {
        std::set<std::string> strSet(buildOrder.begin(), buildOrder.end());
        return strSet.size();
    }


public:

    typedef T race1;
    typedef U race2;

    MOGA(Vec2D const minPos, Vec2D const maxPos, std::string const& filePath1, std::string const& filePath2, size_t popSize, std::vector<std::string> const & buildList1, std::vector<std::string> const & buildList2, size_t const nGoals)
        :  mPopSize(popSize), mDistribution(0,1.0), mNGoals(nGoals)
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
        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
        generator.discard(1000);

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

        mMutationProbability = 1.0/mNGenes;

        if(mNCrossoverPoints > 1) --mNCrossoverPoints;
        std::uniform_real_distribution<double> valueDist(MIN, MAX);
        std::vector<Chromosome> initChroms(mNGoals);
        for(auto& chrom : initChroms)
        {
            chrom.resize(mNGenes);
            for(auto& gene : chrom)
            {
                gene = valueDist(generator);
            }
        }
        setGoals(initChroms);


        mPop.reserve(2*popSize);
        mPop.resize(popSize);


        mStats.first.sum = 0.0;
        mStats.first.max = 0.0;
        mStats.second.sum = 0.0;
        mStats.second.max = 0.0;

        for(size_t i = 0; i < popSize; ++i)
        {
            Chromosome& chrom = mPop[i].chromosome;
            chrom.resize(mNGenes);

            for(auto& gene : chrom)
            {
                gene = valueDist(generator);
            }

        }
        evaluate(mPop);
        mStats.first.onlinePerformance = 0.0;
        mStats.first.offlinePerformance = 0.0;
        mStats.second.onlinePerformance = 0.0;
        mStats.second.offlinePerformance = 0.0;

    }

    void optimize(std::vector<Chromosome> const& goals, size_t const iterations, std::mt19937& generator)
    {
        mMutationCount = 0;
        mIndividualToMutate = 0;
        mGeneToMutate = 0;

        setGoals(goals);
        evaluate(mPop);

        std::vector<double> v(mPop.size());
        for(size_t i = 0; i < v.size(); ++i)
        {
            v[i] = mPop[i].fitness.damage;
        }
        computeStatistics(mStats.first, v);
        for(size_t i = 0; i < v.size(); ++i)
        {
            v[i] = mPop[i].fitness.health;
        }
        computeStatistics(mStats.second, v);

        for(size_t i = 0; i < iterations; ++i)
        {

            mutationClock(generator);
            // apply genetic algorithm

            std::vector<Individual> newPop;
            newPop.reserve(mPopSize);
            std::vector<Individual *> selected;

            selected = tournamentSelection(mPopSize, generator);
            if(selected.size() == 0) break; // Just for safety
            #pragma omp parallel
            {

                std::vector<Individual> newPop_local;
                newPop_local.reserve(selected.size() / 2 / omp_get_num_threads());
                #pragma omp for schedule(dynamic,1) nowait
                for(size_t j = 0; j < selected.size()-1; j += 2)
                {
                    std::pair<Individual, Individual> children;
                    children = crossoverFuncs[mCrossoverChoice](CrossoverParameter({*selected[j], *selected[j+1]}, generator, mNCrossoverPoints));
                    newPop_local.push_back(children.first);
                    newPop_local.push_back(children.second);
                }
                #pragma omp critical
                {
                    newPop.insert(newPop.end(), newPop_local.begin(), newPop_local.end());
                }
            }


            mMutationCount += newPop.size();
            size_t pos = 0, pos_old;
            while(mIndividualToMutate < mMutationCount)
            {
                pos_old = pos;
                size_t const offset = mMutationCount - mIndividualToMutate;
                pos = newPop.size() > offset ? newPop.size() - offset : 0;
                mutationFuncs[mMutationChoice](MutationParameter(newPop[pos], generator, mGeneToMutate, i, iterations));
                mutationClock(generator);
                size_t const pos_diff = pos - pos_old;
                mMutationCount = mMutationCount > pos_diff ? mMutationCount - pos_diff : 0;
            }
            evaluate(newPop);
            mPop.insert(mPop.begin(), newPop.begin(), newPop.end());

            nondominatedSort();
            while(mPop.size() > mPopSize)
            {
                mPop.pop_back();
            }

            std::vector<double> v(mPop.size());
            for(size_t i = 0; i < v.size(); ++i)
            {
                v[i] = mPop[i].fitness.damage;
            }
            computeStatistics(mStats.first, v);
            for(size_t i = 0; i < v.size(); ++i)
            {
                v[i] = mPop[i].fitness.health;
            }
            computeStatistics(mStats.second, v);
        }


    }


    size_t getNumberOfCrossoverOperators() const
    {
        return crossoverFuncs.size();
    }

    size_t getNumberOfMutationOperators() const
    {
        return mutationFuncs.size();
    }

    std::pair<Statistics,Statistics> getStatistics() const
    {
        return mStats;
    }


    void setCrossover(size_t const value)
    {
        mCrossoverChoice = value < crossoverFuncs.size() ? value : mCrossoverChoice;
    }

    void setMutation(size_t const value)
    {
        mMutationChoice = value < mutationFuncs.size() ? value : mMutationChoice;
    }

    std::string getCrossoverOperatorName()
    {
        return mCrossoverFuncNames[mCrossoverChoice];
    }
    std::string getMutationOperatorName()
    {
        return mMutationFuncNames[mMutationChoice];
    }

    std::vector<Chromosome> getBestChromosomes(size_t const n)
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
        std::vector<Chromosome> res;
        res.reserve(sz);
        for(size_t i = 0; i < sz; ++i)
        {
            res.push_back(mPop[i].chromosome);
        }
        return res;
    }

    std::vector<Individual> const& getPopulation() const
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
    size_t getNumberOfGoals()
    {
        return mNGoals;
    }

    void setNumberOfGoals(size_t const value)
    {
        mNGoals = value;
    }

    void includeDecodedChromosomes(Chromosome const& data, size_t const migrants, int const rank, int const procs)
    {
        std::vector<Individual> newPop;
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
            mPop.pop_back();
        } while(mPop.size() > mPopSize);
    }

    unsigned long getNumberOfEvaluations() const
    {
        unsigned long sum = 0;
        for(auto const& simsPerThread : mSims)
        {
            for(auto const& sim : simsPerThread)
            {
                sum += sim.getNumberOfRuns();
            }
        }
        return sum;
    }

    void setTournamentSize(size_t const value)
    {
        mTournamentSize = value > 1 && value < mPopSize/2 ? value : 4;
    }

    size_t getTournamentSize() const
    {
        return mTournamentSize;
    }




};


#endif // MOGA

