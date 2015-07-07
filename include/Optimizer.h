#ifndef _OPTIMIZER_
#define _OPTIMIZER_

#include <ostream>
#include <iostream>
#include <mpi.h>


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
    size_t NGoals;

    MicroSimulation<typename GA1::race1, typename GA2::race1> mSim;


public:
    Optimizer(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2, size_t popSize, vector<string> const & buildList1, vector<string> const & buildList2, size_t const nGoals)
        : popSize(popSize),
          ga1(minPos, maxPos, filePath1, filePath2, popSize, buildList1, buildList2, nGoals),
          ga2(minPos, maxPos, filePath1, filePath2, popSize, buildList1, buildList2, nGoals),
          NGoals(nGoals),
          mSim(minPos, maxPos, filePath1, filePath2)
    {
        mSim.initBothPlayers(buildList1, buildList2);
    }

    void optimize(size_t const sel, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt, int const rank, int const procs, size_t const migrants)
    {
        setSelection(sel);
        setCrossover(co);
        setMutation(mut);

        if(rank == 0)
        {
            cout << "Performing Optimization with the following parameters" << endl << endl;
            cout << "Selection Operator: " << getSelectionOperatorName() << endl;
            cout << "Crossover Operator: " << getCrossoverOperatorName() << endl;
            cout << "Mutation Operator: " << getMutationOperatorName() << endl;
            cout << "Population Size: " << popSize << endl;
            cout << "Number of Iterations: " << iterations << endl;
            cout << "Generations per Iteration: " << genPerIt << endl;
        }

        vector<unsigned long> buf1(procs*migrants*ga1.getNumberOfGenes());
        vector<unsigned long> buf2(procs*migrants*ga2.getNumberOfGenes());

        for(size_t i = 0; i < iterations; ++i)
        {
//            std::cout << "Progress: " << static_cast<double>(i)/iterations*100 << "%" << "\r" << std::flush;
//            printf("%c[2K", 27);
            vector<Chromosome> optima1(ga1.getBestChromosomes(NGoals));
            vector<Chromosome> optima2(ga2.getBestChromosomes(NGoals));
            ga1.optimize(optima2, genPerIt);
            ga2.optimize(optima1, genPerIt);
            stats1 = ga1.getStatistics();
            stats2 = ga2.getStatistics();
            onlinePerformance += 0.5*(ga1.getOnlinePerformance() + ga2.getOnlinePerformance());
            offlinePerformance += 0.5*(ga1.getOfflinePerformance() + ga2.getOfflinePerformance());

            migrate(migrants, buf1, ga1, rank);
            migrate(migrants, buf2, ga2, rank);
        }
        onlinePerformance /= iterations;
        offlinePerformance /= iterations;

        double onlinePerformance_tmp, offlinePerformance_tmp;

        MPI_Reduce (&onlinePerformance, &onlinePerformance_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce (&offlinePerformance, &offlinePerformance_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        computeGlobalStatistics(stats1, rank, procs);
        computeGlobalStatistics(stats2, rank, procs);

        if(rank == 0)
        {
            onlinePerformance = onlinePerformance_tmp/procs;
            offlinePerformance = offlinePerformance_tmp/procs;

            printStatistics();
            cout << "\n" << endl;
        }
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

    bool determineWinner(std::ostream& stream, int const rank, int const procs)
    {
        gatherPopulation(ga1, rank, procs);
        gatherPopulation(ga2, rank, procs);
        if(rank == 0)
        {
            vector<Individual> pop1(ga1.getPopulation());
            vector<Individual> pop2(ga2.getPopulation());
            size_t const minSize = std::min(pop1.size(), pop2.size());
            pop1.resize(minSize);
            pop2.resize(minSize);
            Fitness res;


            for(size_t i = 0; i < minSize; ++i)
            {
                for(size_t j = 0; j < minSize; ++j)
                {
                    mSim.setPlayer1Chromosome(pop1[i].chromosome);
                    mSim.setPlayer2Chromosome(pop2[j].chromosome);
                    res += mSim.run(true, Player::first);
                }
            }
            double const damage1 = (res.damage/(minSize * minSize));
            double const damage2 = (1.0-(res.health/(minSize * minSize)));
            stream << "Comparison of the final populations" << std::endl;
            stream << "Damage caused by Player 1: " << damage1*100 << " %" << std::endl;
            stream << "Remaining health for Player 1: " << (1.0-damage2)*100 << " %" << std::endl;
            stream << "Damage caused by Player 2: " << damage2*100 << " %" << std::endl;
            stream << "Remaining health of Player 2: " << (1.0-damage1)*100 << " %" << std::endl;
            return damage1 > damage2;
        }
    }

    template<typename GA>
    void migrate(size_t const n, vector<unsigned long>& buf, GA& ga, int const rank)
    {
        vector<unsigned long> sendData = std::move(ga.getDecodedChromosomes(n));
        MPI_Allgather(sendData.data(), sendData.size(), MPI_UNSIGNED_LONG, buf.data(), sendData.size(), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        ga.includeDecodedChromosomes(buf, rank);
    }

    void computeGlobalStatistics(Statistics& stats, int const rank, int const procs)
    {
        Statistics tmp;
        MPI_Reduce (&stats.mean, &tmp.mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce (&stats.max, &tmp.max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce (&stats.sum, &tmp.sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce (&stats.stdev, &tmp.stdev, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if(rank == 0)
        {
            stats.mean = tmp.mean / procs;
            stats.max = tmp.max;
            stats.sum = tmp.sum;
            stats.stdev = tmp.stdev / procs;
        }
    }

    template<typename GA>
    void gatherPopulation(GA& ga, int const rank, int const procs)
    {
        vector<unsigned long> buf;
        if(rank == 0)
        {
            buf.resize(procs*popSize*ga.getNumberOfGenes());
        }
        vector<unsigned long> sendData = std::move(ga.getDecodedChromosomes(popSize));
        MPI_Gather(sendData.data(), sendData.size(), MPI_UNSIGNED_LONG, buf.data(), sendData.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        if(rank == 0)
        {
            ga.includeDecodedChromosomes(buf, 0);
        }
    }
};

#endif
