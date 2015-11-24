#include "../include/OptimizerInterface.h"


template<typename Race1, typename Race2>
OptimizerInterface<Race1,Race2>::OptimizerInterface(Vec2D const& minPos, Vec2D const& maxPos, std::string const& filePath1, std::string const& filePath2, size_t popSize, std::vector<std::string> const & buildOrder1, std::vector<std::string> const & buildOrder2, size_t const nGoals, std::string const& resDirPath )
    : mOpt(minPos, maxPos, filePath1, filePath2, popSize, buildOrder1, buildOrder2, nGoals, resDirPath) {}

template<typename Race1, typename Race2>
void OptimizerInterface<Race1,Race2>::optimize(size_t const sel, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt, int const rank, int const procs, size_t migrants, bool saveStatistics)
{
    mOpt.optimize(sel, co, mut, iterations, genPerIt, rank, procs, migrants, saveStatistics);
}

template<typename Race1, typename Race2>
bool OptimizerInterface<Race1,Race2>::determineWinner(std::ostream& stream, int const rank, int const procs, bool saveStatistics)
{
    mOpt.determineWinner(stream, rank, procs, saveStatistics);
}

