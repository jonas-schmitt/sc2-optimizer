#ifndef CHROMOSOME
#define CHROMOSOME
#include<bitset>
#include<vector>
#include<cmath>


typedef std::vector<double> Chromosome;


struct Fitness final
{

    double score = 0.0;
    double damage = 0.0;
    double health = 0.0;
    double timeSteps = 0.0;


    Fitness& operator+=(Fitness const& rhs)
    {
        score += rhs.score;
        damage += rhs.damage;
        health += rhs.health;
        timeSteps += rhs.timeSteps;

        return *this;
    }

    Fitness& operator=(double const value)
    {
        score = value;
        damage = value;
        health = value;
        timeSteps = value;
        return *this;
    }

    Fitness& operator*=(double const value)
    {
        score *= value;
        damage *= value;
        health *= value;
        timeSteps *= value;

        return *this;
    }
    Fitness& operator/=(double const value)
    {
        score /= value;
        damage /= value;
        health /= value;
        timeSteps /= value;
        return *this;
    }

    bool dominates(Fitness const& other) const
    {
        if(this->damage > other.damage && this->health > other.health)
        {
            return true;
        }
        else
        {
            return false;
        }
    }


};

// Struct containing all information associated to an individual
struct Individual
{
    Fitness fitness;
    Chromosome chromosome;
    Chromosome alternative;
    double cdf;
    double total;

    int dominationCount;
    std::vector<size_t> dominationSet;
    int rank;
    double distance;
    bool evaluated = false;
    bool extinction = false;


    Individual() {}
    Individual(size_t const N) : chromosome(N) {}

    bool dominates(Individual const& ind) const
    {
        return fitness.dominates(ind.fitness);
    }

};


#endif // CHROMOSOME

