#ifndef _OPTIMIZER_
#define _OPTIMIZER_

#include <ostream>
#include <iostream>
#include <mpi.h>


#include "moga.h"
#include "Chromosome.h"

// Class implementing the optimization algorithm
// Not intended for direct use

template <typename GA1, typename GA2>
class Optimizer
{
private:
    size_t mPopSize;

    // Genetic algorithm objects
    GA1 mGa1;
    GA2 mGa2;

    // Statistics of both player's populations
    std::pair<Statistics,Statistics> mStats1, mStats2;
    double mOfflinePerformance = 0.0;
    double mOnlinePerformance = 0.0;

    // Number of strategies used for fitness evaluation
    size_t mNGoals;

    // Multiple colonies?
    bool mMPI;

    //
    std::string mResDirPath;


    // Simulation used for the final comparison
    std::vector<MicroSimulation<typename GA1::race1, typename GA2::race2>> mSim;


public:
    Optimizer(Vec2D const& minPos, Vec2D const& maxPos, std::string const& filePath1, std::string const& filePath2, size_t popSize, std::vector<std::string> const & buildOrder1, std::vector<std::string> const & buildOrder2, size_t const nGoals, std::string const& resDirPath)
        : mPopSize(popSize),
          mGa1(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals),
          mGa2(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals),
          mNGoals(nGoals),
          mResDirPath(resDirPath)
    {
        size_t nthreads = 1;
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                nthreads = omp_get_num_threads();
            }
        }
        for(size_t i = 0; i < nthreads; ++i)
        {
            mSim.emplace_back(minPos, maxPos, filePath1, filePath2);
        }
        for(size_t i = 0; i < nthreads; ++i)
        {
            mSim[i].initBothPlayers(buildOrder1, buildOrder2);
        }
    }

    void optimize(size_t const tournamentSize, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt, int const rank, int const procs, size_t migrants, bool saveStatistics)
    {
        static std::mt19937 generator;

        #pragma omp threadprivate(generator)
        #pragma omp parallel
        {
            generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
            generator.discard(1000);
        }

        mMPI = procs > 1;
        setTournamentSize(tournamentSize);
        setCrossover(co);
        setMutation(mut);

        migrants = std::min(migrants, mPopSize);

        if(rank == 0)
        {
            std::cout << "Performing Optimization with the following parameters" << std::endl << std::endl;
            std::cout << "Tournament Selection with a Tournament Size of: " << getTournamentSize() << std::endl;
            std::cout << "Crossover Operator: " << getCrossoverOperatorName() << std::endl;
            std::cout << "Mutation Operator: " << getMutationOperatorName() << std::endl;
            std::cout << "Population Size: " << mPopSize*procs << std::endl;
            std::cout << "Number of Iterations: " << iterations << std::endl;
            std::cout << "Generations per Iteration: " << genPerIt << std::endl;
        }

        Chromosome buf1(procs*migrants*mGa1.getNumberOfGenes());
        Chromosome buf2(procs*migrants*mGa2.getNumberOfGenes());
        std::ofstream avgFile1;
        std::ofstream avgFile2;
        std::ofstream stdevFile1;
        std::ofstream stdevFile2;
        if(saveStatistics)
        {
            avgFile1.open(mResDirPath+"/avgPlayer1.dat");
            avgFile2.open(mResDirPath+"/avgPlayer2.dat");
            stdevFile1.open(mResDirPath+"/stdevPlayer1.dat");
            stdevFile2.open(mResDirPath+"/stdevPlayer2.dat");
            mGa1.writeOutStatistics(avgFile1, stdevFile1);
            mGa2.writeOutStatistics(avgFile2, stdevFile2);
        }

        mGa1.setNumberOfSamples(2*genPerIt+1);
        mGa2.setNumberOfSamples(2*genPerIt+1);
        // Perform the optimization competitively for both players for iterations x genPerIt generations
        for(size_t i = 0; i < iterations; ++i)
        {

            if(rank == 0)
            {
                std::cout << "Progress: " << static_cast<double>(i)/iterations*100 << " %" << "\r" << std::flush;
                printf("%c[2K", 27);
            }

            // Migration
            if(mMPI)
            {
                migrate(buf1, migrants, mGa1, rank, procs);
                migrate(buf2, migrants, mGa2, rank, procs);
            }

            std::vector<Chromosome> optima1(mGa1.getBestChromosomes(mNGoals));
            std::vector<Chromosome> optima2(mGa2.getBestChromosomes(mNGoals));

            // Optimize using the mNGoals best chromosomes of the opponent for fitness evaluation
            bool const retVal1 = mGa1.optimize(optima2, genPerIt, generator, rank, procs);
            bool const retVal2 = mGa2.optimize(optima1, genPerIt, generator, rank, procs);

            // Get process-local statistics
            mStats1 = mGa1.getStatistics();
            mStats2 = mGa2.getStatistics();

            if(mMPI)
            {
                computeGlobalStatistics(mStats1.first, rank, procs);
                computeGlobalStatistics(mStats1.second, rank, procs);
                computeGlobalStatistics(mStats2.first, rank, procs);
                computeGlobalStatistics(mStats2.second, rank, procs);
            }
            if(retVal1 && retVal2)
            {
                break;
            }

        }
        unsigned long evaluations = mGa1.getNumberOfEvaluations() + mGa2.getNumberOfEvaluations();
        unsigned long evaluations_tmp;
        if(mMPI)
        {
            MPI_Reduce (&evaluations, &evaluations_tmp, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        if(rank == 0)
        {
            if(mMPI)
            {
                evaluations = evaluations_tmp;
            }

            std::cout << "Number of evaluations: " << evaluations << std::endl;
            std::cout << std::endl;
            std::cout << "Damage" << std::endl;
            std::cout << "Overall Online Performance: " << (mStats1.first.onlinePerformance + mStats2.first.onlinePerformance)/(mStats1.first.iteration + mStats2.first.iteration)*100.0 << std::endl;
            std::cout << "Overall Offline Performance: " << (mStats1.first.offlinePerformance + mStats2.first.offlinePerformance)/(mStats1.first.iteration + mStats2.first.iteration)*100.0 << std::endl;
            std::cout << std::endl;
            std::cout << "Health" << std::endl;
            std::cout << "Overall Online Performance: " << (mStats1.second.onlinePerformance + mStats2.second.onlinePerformance)/(mStats1.second.iteration + mStats2.second.iteration)*100.0 << std::endl;
            std::cout << "Overall Offline Performance: " << (mStats1.second.offlinePerformance + mStats2.second.offlinePerformance)/(mStats1.second.iteration + mStats2.second.iteration)*100.0 << std::endl;
            std::cout << "\n" << std::endl;
        }
        if(saveStatistics)
        {
            mGa1.stopWriteOutStatistics();
            mGa2.stopWriteOutStatistics();
            avgFile1.close();
            avgFile2.close();
            stdevFile1.close();
            stdevFile2.close();
        }
    }

    void printStatistics()
    {
        std::string separator;
        for(int i = 0; i < 50; ++i) separator += "-";
        separator += '\n';
        std::cout << separator << "Statistics Player 1:" << std::endl;
        std::cout << "Damage" << std::endl;
        mStats1.first.print();
        std::cout << "Health" << std::endl;
        mStats1.second.print();
        std::cout << separator << std::endl;
        std::cout << separator << "Statistics Player 2:" << std::endl;
        std::cout << "Damage" << std::endl;
        mStats2.first.print();
        std::cout << "Health" << std::endl;
        mStats2.second.print();
        std::cout << separator << std::endl;
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

    void setTournamentSize(size_t const value)
    {
        mGa1.setTournamentSize(value);
        mGa2.setTournamentSize(value);
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

    size_t getTournamentSize()
    {
        return mGa1.getTournamentSize();
    }
    std::string getCrossoverOperatorName()
    {
        return mGa1.getCrossoverOperatorName();
    }
    std::string getMutationOperatorName()
    {
        return mGa1.getMutationOperatorName();
    }

    bool determineWinner(std::ostream& stream, int const rank, int const procs, bool saveStatistics)
    {
        // Gather the combined population of all colonies
        if(mMPI)
        {
            gatherPopulation(mGa1, rank, procs);
            gatherPopulation(mGa2, rank, procs);
        }

        if(rank == 0)
        {

            // Final comparison of the 100 best individuals of the opponents
            std::vector<Individual> pop1(mGa1.getPopulation());
            std::vector<Individual> pop2(mGa2.getPopulation());
            auto cmp = [] (Individual const& lhs, Individual const& rhs)
            {
                return lhs.fitness.score > rhs.fitness.score;
            };
//            std::stable_sort(pop1.begin(), pop1.end(), cmp);
//            std::stable_sort(pop2.begin(), pop2.end(), cmp);



            size_t const samples = std::min(static_cast<size_t>(1000UL), std::min(pop1.size(), pop2.size()));

            pop1.resize(samples);
            pop2.resize(samples);
            for(size_t i = 0; i < samples; ++i)
            {
                pop1[i].fitness = 0;
                pop2[i].fitness = 0;
            }

            Fitness res;

            // Increase the time step size
            std::cout << "Generating final result" << std::endl;
            size_t progress = 0;
            #pragma omp parallel for schedule(static)
            for(size_t i = 0; i < samples; ++i)
            {
                for(size_t j = 0; j < samples; ++j)
                {
                    mSim[omp_get_thread_num()].setPlayer1Chromosome(pop1[i].chromosome);
                    mSim[omp_get_thread_num()].setPlayer2Chromosome(pop2[j].chromosome);
                    Fitness tmp1 = mSim[omp_get_thread_num()].run(true, Player::first);
                    Fitness tmp2;
                    tmp2.damage = 1.0 - tmp1.health;
                    tmp2.health = 1.0 - tmp1.damage;
                    tmp2.score = 50.0*(tmp2.damage + tmp2.health);
                    tmp1 /= samples;
                    pop1[i].fitness += tmp1;
                    tmp2 /= samples;
                    #pragma omp critical
                    {
                        res += tmp1;
                        pop2[j].fitness += tmp2;
                    }
                }
                #pragma omp critical
                {
                    ++progress;
                    std::cout << "Progress: " << static_cast<double>(progress)/static_cast<double>(samples)*100.0 << " %" << "\r" << std::flush;
                    printf("%c[2K", 27);
                }
            }

            // Sort the individuals of the final comparison according to their performance
            std::sort(pop1.begin(), pop1.end(), cmp);
            std::sort(pop2.begin(), pop2.end(), cmp);

            // Compute and store the result
            double const damage1 = (res.damage/(samples));
            double const damage2 = (1.0-(res.health/(samples)));
//            stream << "Average of Player 1: " << "\n";
//            stream << "Damage: " << mStats1.first.mean << "\n";
//            stream << "Health: " << mStats1.second.mean << "\n";
//            stream << "Average of Player 2: " << "\n";
//            stream << "Damage: " << mStats2.first.mean << "\n";
//            stream << "Health: " << mStats2.second.mean << "\n";
            stream << "Comparison of the final populations" << std::endl;
            stream << "Damage caused by Player 1: " << damage1*100 << " %" << std::endl;
            stream << "Remaining health for Player 1: " << (1.0-damage2)*100 << " %" << std::endl;
            stream << "Damage caused by Player 2: " << damage2*100 << " %" << std::endl;
            stream << "Remaining health of Player 2: " << (1.0-damage1)*100 << " %" << std::endl << std::endl;



            if(saveStatistics)
            {
                for(size_t i = 0; i < std::min(static_cast<size_t>(10), samples); ++i)
                {
                    for(size_t j = 0; j < std::min(static_cast<size_t>(10), samples); ++j)
                    {
                        mSim.front().setPlayer1Chromosome(pop1[i].chromosome);
                        mSim.front().setPlayer2Chromosome(pop2[j].chromosome);
                        mSim.front().enableTracking(mResDirPath+"/pl1_paths_" + std::to_string(i) + "_" + std::to_string(j) + ".dat",
                                            mResDirPath+"/pl2_paths_" + std::to_string(i) + "_" + std::to_string(j) + ".dat");
                        mSim.front().run(true, Player::first);
                        mSim.front().disableTracking();
                    }
                }
            }
        }
    }

    // Migration of migrants individuals
    // buf: Buffer for receiving the data
    template<typename GA>
    void migrate(Chromosome& buf, size_t const migrants, GA& ga, int const rank, int const procs)
    {
        // Contains the concatenated chromosomes of multiple individuals
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
        MPI_Reduce (&stats.onlinePerformance, &tmp.onlinePerformance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce (&stats.offlinePerformance, &tmp.offlinePerformance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if(rank == 0)
        {
            stats.mean = tmp.mean / procs;
            stats.max = tmp.max;
            stats.sum = tmp.sum;
            stats.stdev = tmp.stdev / procs;
            stats.onlinePerformance = tmp.onlinePerformance / procs;
            stats.offlinePerformance = tmp.offlinePerformance / procs;
        }
    }

    // Gather the population of all colonies at the root process
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
