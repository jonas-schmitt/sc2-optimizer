#ifndef OPTIMIZERINTERFACE
#define OPTIMIZERINTERFACE

#include <string>
#include <vector>
#include <fstream>
#include "Optimizer.h"
#include "moga.h"
#include "Race.h"

template<typename Race1, typename Race2> class OptimizerInterface
{
private:
    Optimizer<MOGA<Race1, Race2, Player::first>, MOGA<Race1, Race2, Player::second> > mOpt;

public:
    OptimizerInterface(Vec2D const& minPos, Vec2D const& maxPos, std::string const& filePath1, std::string const& filePath2, size_t popSize, std::vector<std::string> const & buildList1, std::vector<std::string> const & buildList2, size_t const nGoals, std::string const& dirPath )
        : mOpt(minPos, maxPos, filePath1, filePath2, popSize, buildList1, buildList2, nGoals, dirPath) {}

    void optimize(size_t const sel, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt, int const rank, int const procs, size_t migrants, bool saveStatistics)
    {
        mOpt.optimize(sel, co, mut, iterations, genPerIt, rank, procs, migrants, saveStatistics);
    }

    bool determineWinner(std::ostream& stream, int const rank, int const procs, bool saveStatistics)
    {
        mOpt.determineWinner(stream, rank, procs, saveStatistics);
    }

};

#endif // OPTIMIZERINTERFACE

