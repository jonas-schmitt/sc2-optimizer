#ifndef CHROMOSOME
#define CHROMOSOME
#include<bitset>
#include<vector>
#include<cmath>


//static size_t const NBITS = 32;

//typedef std::vector<bitset<NBITS>> Chromosome;
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


    Individual() {}
    Individual(size_t const N) : chromosome(N) {}

    size_t computeHash() const
    {
        size_t hash = 0;
        //std::hash<bitset<NBITS> > hash_func;
        std::hash<double > hash_func;
        for(auto const gene : chromosome)
        {
            hash ^= hash_func(gene) + 0x9e3779b9 + (hash<<6) + (hash>>2);
        }
        return hash;
    }

    bool dominates(Individual const& ind) const
    {
        return fitness.dominates(ind.fitness);
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

