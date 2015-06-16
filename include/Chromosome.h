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
    double damage_killed = 0.0;
    double gas_killed = 0.0;
    double minerals_killed = 0.0;

    double health = 0.0;
    double health_alive = 0.0;
    double gas_alive = 0.0;
    double minerals_alive = 0.0;



};

struct Individual
{
    Fitness fitness;
    Chromosome chromosome;
    double cdf;
    Individual() {}
    Individual(size_t const N) : chromosome(N) {}
};




#endif // CHROMOSOME

