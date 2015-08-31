#ifndef OPTIMIZERINTERFACE
#define OPTIMIZERINTERFACE
#include "Optimizer_real.h"
#include "moga_real.h"
#include "Race.h"

template<typename Race1, typename Race2> class OptimizerInterface
{
private:
    Optimizer_real<MOGA_real<Race1, Race2, Player::first>, MOGA_real<Race1, Race2, Player::second> > mOpt;

public:
    OptimizerInterface(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2, size_t popSize, vector<string> const & buildList1, vector<string> const & buildList2, size_t const nGoals)
        : mOpt(minPos, maxPos, filePath1, filePath2, popSize, buildList1, buildList2, nGoals) {}

    void optimize(size_t const sel, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt, int const rank, int const procs, size_t migrants)
    {
        mOpt.optimize(sel, co, mut, iterations, genPerIt, rank, procs, migrants);
    }

    bool determineWinner(std::ostream& stream, int const rank, int const procs)
    {
        mOpt.determineWinner(stream, rank, procs);
    }

};

#endif // OPTIMIZERINTERFACE

