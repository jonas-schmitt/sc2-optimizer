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
    double health = 0.0;


    Fitness& operator+=(Fitness const& rhs)
    {
        score += rhs.score;
        damage += rhs.damage;
        health += rhs.health;

        return *this;
    }

    Fitness& operator=(double const value)
    {
        score = value;
        damage = value;
        health = value;
        return *this;
    }

    Fitness& operator*=(double const value)
    {
        score *= value;
        damage *= value;
        health *= value;

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
        if(ind.fitness.damage - fitness.damage < EPS
                && ind.fitness.health - fitness.health < EPS)
        {
            return true;
        }
        else
        {
            return false;
        }
    }


};

//inline bool operator< (Individual const& lhs, Individual const& rhs)
//{
//    if(lhs.rank < rhs.rank)
//    {
//        return true;
//    }
//    else if(lhs.rank == rhs.rank)
//    {
//        return lhs.distance > rhs.distance;
//    }
//    else
//    {
//        return false;
//    }
//}




#endif // CHROMOSOME

