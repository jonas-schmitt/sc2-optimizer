#ifndef OPTIMIZERINTERFACE
#define OPTIMIZERINTERFACE

#include <string>
#include <vector>
#include <fstream>
#include "Optimizer.h"
#include "moga.h"
#include "Race.h"


// Interface for using the Optimizer class

template<typename Race1, typename Race2> class OptimizerInterface
{
private:
    Optimizer<MOGA<Race1, Race2, Player::first>, MOGA<Race1, Race2, Player::second> > mOpt;

public:

    // Initialize the optimizer
    // minPos, maxPos: Minimal and maximal coordinates of the playground
    // filePath1, filePath2: Paths to the directories containing the unit stats
    // popSize: Population size
    // buildOrder1, buildOrder2: Build orders of both players, denoted by the names of the units produced
    // nGoals: Number of strategies used for fitness evaluation
    // resDirPath: Path to the directory where the result files should be saved
    OptimizerInterface(Vec2D const& minPos, Vec2D const& maxPos, std::string const& filePath1, std::string const& filePath2, size_t popSize, std::vector<std::string> const & buildOrder1, std::vector<std::string> const & buildOrder2, size_t const nGoals, std::string const& resDirPath);

    // Run the optimization
    // sel: Tournament size for selection
    // co: Crossover operator choice
    // mut: Mutation operator choice
    // iterations x genPerIt: Number of total generations
    // rank: Rank of the process
    // procs: Total number of processes/colonies
    // migrants: Number of migrants
    // saveStatistics: Should additional statistical information about the optimization be saved?
    void optimize(size_t const sel, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt, int const rank, int const procs, size_t migrants, bool saveStatistics);

    // Determine the winner of the encounter
    // stream: File stream for outputting the result
    // rank, procs, saveStatistics: As above
    bool determineWinner(std::ostream& stream, int const rank, int const procs, bool saveStatistics);

};

#endif // OPTIMIZERINTERFACE

