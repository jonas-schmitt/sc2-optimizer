#ifndef CHROMOSOME
#define CHROMOSOME
#include<bitset>
#include<vector>
#include<cmath>


using std::vector;
using std::bitset;

static size_t const NBITS = 32;

typedef vector<bitset<NBITS>> Chromosome;


struct Fitness final
{

    double score = 0.0;
    double damage = 0.0;
    double gas_killed = 0.0;
    double minerals_killed = 0.0;

    double health = 0.0;
    double gas_alive = 0.0;
    double minerals_alive = 0.0;



};

struct Individual
{
    Fitness fitness;
    Chromosome chromosome;
    double cdf;
    double total;
    Individual() {}
    Individual(size_t const N) : chromosome(N) {}

    size_t computeHash() const
    {
        size_t hash = 0;
        std::hash<bitset<NBITS> > hash_func;
        for(auto const bits : chromosome)
        {
            hash ^= hash_func(bits) + 0x9e3779b9 + (hash<<6) + (hash>>2);
        }
        return hash;
    }
};




#endif // CHROMOSOME

