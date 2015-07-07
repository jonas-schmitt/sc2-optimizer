#ifndef _OPTIMIZER_
#define _OPTIMIZER_

#include <ostream>
#include <iostream>


#include "soga.h"
#include "moga.h"
#include "Chromosome.h"

using std::cout;
using std::endl;

template <typename GA1, typename GA2>
class Optimizer
{
private:
    size_t popSize;
    GA1 ga1;
    GA2 ga2;
    Individual optimum1, optimum2;
    Statistics stats1, stats2;
    double offlinePerformance = 0.0;
    double onlinePerformance = 0.0;

    MicroSimulation<typename GA1::race1, typename GA2::race1> mSim;


public:
    Optimizer(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2, size_t popSize, vector<string> const & buildList1, vector<string> const & buildList2, size_t const nGoals)
        : popSize(popSize),
          ga1(minPos, maxPos, filePath1, filePath2, popSize, buildList1, buildList2, nGoals),
          ga2(minPos, maxPos, filePath2, filePath1, popSize, buildList2, buildList1, nGoals),
          mSim(minPos, maxPos, filePath1, filePath2)
    {
        mSim.initBothPlayers(buildList1, buildList2);
    }

    void optimize(size_t const sel, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt)
    {
        setSelection(sel);
        setCrossover(co);
        setMutation(mut);

        cout << "Performing Optimization with the following parameters" << endl << endl;
        cout << "Selection Operator: " << getSelectionOperatorName() << endl;
        cout << "Crossover Operator: " << getCrossoverOperatorName() << endl;
        cout << "Mutation Operator: " << getMutationOperatorName() << endl;
        cout << "Population Size: " << popSize << endl;
        cout << "Number of Iterations: " << iterations << endl;
        cout << "Generations per Iteration: " << genPerIt << endl;

        for(size_t i = 0; i < iterations; ++i)
        {
//            std::cout << "Progress: " << static_cast<double>(i)/iterations*100 << "%" << "\r" << std::flush;
//            printf("%c[2K", 27);
            vector<Chromosome> optima1(ga1.getBestChromosomes());
            vector<Chromosome> optima2(ga2.getBestChromosomes());
            ga1.optimize(optima2, genPerIt);
            ga2.optimize(optima1, genPerIt);
            stats1 = ga1.getStatistics();
            stats2 = ga2.getStatistics();
            onlinePerformance += 0.5*(ga1.getOnlinePerformance() + ga2.getOnlinePerformance());
            offlinePerformance += 0.5*(ga1.getOfflinePerformance() + ga2.getOfflinePerformance());
        }
        onlinePerformance /= iterations;
        offlinePerformance /= iterations;



        printStatistics();
        cout << "\n" << endl;
    }

    void printStatistics()
    {
        string separator;
        for(int i = 0; i < 50; ++i) separator += "-";
        separator += '\n';
        cout << separator << "Statistics Player 1:" << endl;
        stats1.print();
        cout << separator << endl;
        cout << separator << "Statistics Player 2:" << endl;
        stats2.print();
        cout << separator << endl;
        cout << "Online Performance: " << onlinePerformance << endl;
        cout << "Offline Performance: " << offlinePerformance << endl;
        cout << separator << endl;
    }

    GA1 const& getGA1() const
    {
        return ga1;
    }
    GA2 const& getGA2() const
    {
        return ga2;
    }

    size_t getNumberOfSelectionOperators()
    {
        ga1.getNumberOfSelectionOperators();
    }

    size_t getNumberOfCrossoverOperators()
    {
        return ga1.getNumberOfCrossoverOperators();
    }

    size_t getNumberOfMutationOperators()
    {
        return ga1.getNumberOfMutationOperators();
    }



    void setSelection(size_t const value)
    {
        ga1.setSelection(value);
        ga2.setSelection(value);
    }

    void setCrossover(size_t const value)
    {
        ga1.setCrossover(value);
        ga2.setCrossover(value);
    }

    void setMutation(size_t const value)
    {
        ga1.setMutation(value);
        ga2.setMutation(value);
    }

    string getSelectionOperatorName()
    {
        return ga1.getSelectionOperatorName();
    }
    string getCrossoverOperatorName()
    {
        return ga1.getCrossoverOperatorName();
    }
    string getMutationOperatorName()
    {
        return ga1.getMutationOperatorName();
    }

    bool determineWinner(std::ostream& stream)
    {
        vector<Individual> pop1(ga1.getPopulation());
        vector<Individual> pop2(ga2.getPopulation());
        size_t const minSize = std::min(pop1.size(), pop2.size());
        pop1.resize(minSize);
        pop2.resize(minSize);
        Fitness res1;
        Fitness res2;

        for(size_t i = 0; i < minSize; ++i)
        {
            for(size_t j = 0; j < minSize; ++j)
            {
                mSim.setPlayer1Chromosome(pop1[i].chromosome);
                mSim.setPlayer2Chromosome(pop2[j].chromosome);
                res1 += mSim.run(true);
                mSim.setPlayer1Chromosome(pop2[j].chromosome);
                mSim.setPlayer2Chromosome(pop1[i].chromosome);
                res2 += mSim.run(true);
            }
        }
        stream << "Result" << std::endl;
        stream << "Damage caused by Player 1: " << res1.damage*100.0/(minSize * minSize) << " %" << std::endl;
        stream << "Damage caused by Player 2: " << res2.damage*100.0/(minSize * minSize) << " %" << std::endl;
        return res1.damage > res2.damage;
    }
};

#endif
