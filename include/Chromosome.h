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

    double score;
    double damage;
    double damage_killed;
    double gas_killed;
    double minerals_killed;

    double health;
    double health_alive;
    double gas_alive;
    double minerals_alive;

    Fitness()
        : damage(0), health(0), score(0)
    {}
    Fitness(double d, double h, double s)
        : damage(d), health(h), score(s)
    {}
    Fitness(Fitness const& sim)
    {
        damage = sim.damage;
        health = sim.health;
        score = sim.score;
    }
    Fitness operator+(Fitness const& other) const
    {
        return Fitness(this->damage+other.damage,this->health+other.health,this->score+other.score);
    }
    Fitness operator/(double value) const
    {
        return Fitness(this->damage/value,this->health/value,this->score/value);
    }


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

