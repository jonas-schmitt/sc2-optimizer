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
    double minerals_killed = 0.0;
    double gas_killed = 0.0;

    double health = 0.0;
    double minerals_alive = 0.0;
    double gas_alive = 0.0;

    Fitness& operator+=(Fitness const& rhs)
    {
        score += rhs.score;
        damage += rhs.damage;
        minerals_killed += rhs.minerals_killed;
        gas_killed += rhs.gas_killed;
        health += rhs.health;
        minerals_alive += rhs.minerals_alive;
        gas_alive += rhs.gas_alive;
        return *this;
    }

    Fitness& operator=(double const value)
    {
        score = value;
        damage = value;
        minerals_killed = value;
        gas_killed = value;
        health = value;
        minerals_alive = value;
        gas_alive = value;
        return *this;
    }

    Fitness& operator*=(double const value)
    {
        score *= value;
        damage *= value;
        minerals_killed *= value;
        gas_killed *= value;
        health *= value;
        minerals_alive *= value;
        gas_alive *= value;
        return *this;
    }


};

struct Individual
{
    Fitness fitness;
    Chromosome chromosome;
    double cdf;
    double total;

    size_t dominationCount;
    vector<size_t> dominationSet;
    size_t rank;
    double distance;


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

    bool dominates(Individual const& ind) const
    {
        return fitness.damage > ind.fitness.damage && fitness.minerals_killed > ind.fitness.minerals_killed && fitness.gas_killed > ind.fitness.gas_killed
                && fitness.health > ind.fitness.health && fitness.minerals_alive > ind.fitness.minerals_killed && fitness.gas_alive > ind.fitness.gas_alive;
    }


};

inline bool operator< (Individual const& lhs, Individual const& rhs)
{
    if(lhs.rank < rhs.rank)
    {
        return true;
    }
    else if(lhs.rank == rhs.rank)
    {
        return lhs.distance > rhs.distance;
    }
    else
    {
        return false;
    }
}




#endif // CHROMOSOME

