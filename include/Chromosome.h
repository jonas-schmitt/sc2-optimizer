#ifndef CHROMOSOME
#define CHROMOSOME
#include<bitset>
#include<vector>
#include<cmath>

using std::vector;
using std::bitset;

template<size_t N = 32>
class Chromosome
{
public:

    vector<bitset<N>> data;
    Chromosome() {}
    Chromosome(size_t const nGenes) : data(nGenes) {}
    double getPhenotype(size_t const pos) const
    {
        double constexpr max = std::pow(2, N) - 1;
        return static_cast<double>(data[pos].to_ulong())/max;
    }
};

#endif // CHROMOSOME

