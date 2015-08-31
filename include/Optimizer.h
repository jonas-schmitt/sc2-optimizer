#ifndef _OPTIMIZER_
#define _OPTIMIZER_

#include <ostream>
#include <iostream>
#include <mpi.h>


#include "moga.h"
#include "Chromosome.h"

using std::cout;
using std::endl;

template <typename GA1, typename GA2>
class Optimizer
{
private:
    size_t mPopSize;
    GA1 mGa1;
    GA2 mGa2;
    Individual mOptimum1, mOptimum2;
    Statistics mStats1, mStats2;
    double mOfflinePerformance = 0.0;
    double mOnlinePerformance = 0.0;
    size_t mNGoals;
    bool mMPI;


    MicroSimulation<typename GA1::race1, typename GA2::race2> mSim;


public:
    Optimizer(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2, size_t popSize, vector<string> const & buildList1, vector<string> const & buildList2, size_t const nGoals)
        : mPopSize(popSize),
          mGa1(minPos, maxPos, filePath1, filePath2, popSize, buildList1, buildList2, nGoals),
          mGa2(minPos, maxPos, filePath1, filePath2, popSize, buildList1, buildList2, nGoals),
          mNGoals(nGoals),
          mSim(minPos, maxPos, filePath1, filePath2)
    {
        mSim.initBothPlayers(buildList1, buildList2);
    }

    void optimize(size_t const sel, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt, int const rank, int const procs, size_t migrants)
    {
        mMPI = procs > 1;
        setSelection(sel);
        setCrossover(co);
        setMutation(mut);

        migrants = std::min(migrants, mPopSize);

        if(rank == 0)
        {
            cout << "Performing Optimization with the following parameters" << endl << endl;
            cout << "Selection Operator: " << getSelectionOperatorName() << endl;
            cout << "Crossover Operator: " << getCrossoverOperatorName() << endl;
            cout << "Mutation Operator: " << getMutationOperatorName() << endl;
            cout << "Population Size: " << mPopSize*procs << endl;
            cout << "Number of Iterations: " << iterations << endl;
            cout << "Generations per Iteration: " << genPerIt << endl;
        }

        Chromosome buf1(procs*migrants*mGa1.getNumberOfGenes());
        Chromosome buf2(procs*migrants*mGa2.getNumberOfGenes());

        for(size_t i = 0; i < iterations; ++i)
        {
            if(rank == 0)
            {
                std::cout << "Progress: " << static_cast<double>(i)/iterations*100 << "%" << endl;
                //std::cout << "Progress: " << static_cast<double>(i)/iterations*100 << "%" << "\r" << std::flush;
                //printf("%c[2K", 27);
            }

            if(mMPI)
            {
                migrate(buf1, migrants, mGa1, rank, procs);
                migrate(buf2, migrants, mGa2, rank, procs);
            }

            vector<Chromosome> optima1(mGa1.getBestChromosomes(mNGoals));
            vector<Chromosome> optima2(mGa2.getBestChromosomes(mNGoals));

            mGa1.optimize(optima2, genPerIt);
            mGa2.optimize(optima1, genPerIt);

            mStats1 = mGa1.getStatistics();
            mStats2 = mGa2.getStatistics();

            mOnlinePerformance += 0.5*(mGa1.getOnlinePerformance() + mGa2.getOnlinePerformance());
            mOfflinePerformance += 0.5*(mGa1.getOfflinePerformance() + mGa2.getOfflinePerformance());
        }
        mOnlinePerformance /= iterations;
        mOfflinePerformance /= iterations;
        double onlinePerformance_tmp, offlinePerformance_tmp;
        if(mMPI)
        {

            MPI_Reduce (&mOnlinePerformance, &onlinePerformance_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (&mOfflinePerformance, &offlinePerformance_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            computeGlobalStatistics(mStats1, rank, procs);
            computeGlobalStatistics(mStats2, rank, procs);
        }

        if(rank == 0)
        {
            if(mMPI)
            {
                mOnlinePerformance = onlinePerformance_tmp/procs;
                mOfflinePerformance = offlinePerformance_tmp/procs;
            }

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
        mStats1.print();
        cout << separator << endl;
        cout << separator << "Statistics Player 2:" << endl;
        mStats2.print();
        cout << separator << endl;
        cout << "Online Performance: " << mOnlinePerformance << endl;
        cout << "Offline Performance: " << mOfflinePerformance << endl;
        cout << separator << endl;
    }

    GA1 const& getGA1() const
    {
        return mGa1;
    }
    GA2 const& getGA2() const
    {
        return mGa2;
    }

    size_t getNumberOfSelectionOperators()
    {
        mGa1.getNumberOfSelectionOperators();
    }

    size_t getNumberOfCrossoverOperators()
    {
        return mGa1.getNumberOfCrossoverOperators();
    }

    size_t getNumberOfMutationOperators()
    {
        return mGa1.getNumberOfMutationOperators();
    }



    void setSelection(size_t const value)
    {
        mGa1.setSelection(value);
        mGa2.setSelection(value);
    }

    void setCrossover(size_t const value)
    {
        mGa1.setCrossover(value);
        mGa2.setCrossover(value);
    }

    void setMutation(size_t const value)
    {
        mGa1.setMutation(value);
        mGa2.setMutation(value);
    }

    string getSelectionOperatorName()
    {
        return mGa1.getSelectionOperatorName();
    }
    string getCrossoverOperatorName()
    {
        return mGa1.getCrossoverOperatorName();
    }
    string getMutationOperatorName()
    {
        return mGa1.getMutationOperatorName();
    }

    bool determineWinner(std::ostream& stream, int const rank, int const procs)
    {
        if(mMPI)
        {
            gatherPopulation(mGa1, rank, procs);
            gatherPopulation(mGa2, rank, procs);
        }

        if(rank == 0)
        {
            vector<Individual> pop1(mGa1.getPopulation());
            vector<Individual> pop2(mGa2.getPopulation());

            size_t const minSize = std::min(pop1.size(), pop2.size());

            pop1.resize(minSize);
            pop2.resize(minSize);

            auto cmp = [] (Individual const& lhs, Individual const& rhs)
            {
                return lhs.fitness.score > rhs.fitness.score;
            };

            std::sort(pop1.begin(), pop1.end(), cmp);
            std::sort(pop2.begin(), pop2.end(), cmp);

            Fitness res;
            for(size_t i = 0; i < minSize; ++i)
            {
                for(size_t j = 0; j < minSize; ++j)
                {
                    mSim.setPlayer1Chromosome(pop1[i].chromosome);
                    mSim.setPlayer2Chromosome(pop2[j].chromosome);
                    if(i < 5 && j < 5) mSim.enableTracking("./results/pl1_paths_" + std::to_string(i) + "_" + std::to_string(j) + ".txt",
                                                           "./results/pl2_paths_" + std::to_string(i) + "_" + std::to_string(j) + ".txt");
                    res += mSim.run(true, Player::first);
                    if(i < 5 && j < 5) mSim.disableTracking();
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
    void migrate(Chromosome& buf, size_t const migrants, GA& ga, int const rank, int const procs)
    {
        Chromosome sendData = std::move(ga.getChromosomes(migrants));
        MPI_Allgather(sendData.data(), sendData.size(), MPI_DOUBLE, buf.data(), sendData.size(), MPI_DOUBLE, MPI_COMM_WORLD);
        ga.includeDecodedChromosomes(buf, migrants, rank, procs);
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
        Chromosome buf;
        if(rank == 0)
        {
            buf.resize(procs*mPopSize*ga.getNumberOfGenes());
        }
        Chromosome sendData = std::move(ga.getChromosomes(mPopSize));

        if(rank == 0)
        {
            MPI_Gather(sendData.data(), sendData.size(), MPI_DOUBLE, buf.data(), sendData.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            ga.includeDecodedChromosomes(buf, mPopSize, 0, procs);
        }
        else
        {
            MPI_Gather(sendData.data(), sendData.size(), MPI_DOUBLE, NULL, sendData.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
};

#endif
