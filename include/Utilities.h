#ifndef UTILITIES_H
#define UTILITIES_H
#include <cmath>

static double const LIMIT = 10000;
static double const EPS = 1e-9;
static double const GASTOMINERALS = 1.5;
static double const STDLEN = 1./std::sqrt(2);
static long const INTLIMIT = std::numeric_limits<int>::max();
static long const BOUND1 = std::pow(2,8)-1;
static long const BOUND2 = std::pow(2,16)-1;
static int const XMIN = 0;
static int const XMAX = static_cast<int>(std::min(INTLIMIT,BOUND1));
static int const YMIN = 0;
static int const YMAX = static_cast<int>(std::min(INTLIMIT,BOUND2));

struct Fitness
{
    double damage;
    double health;
    double score;
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
#endif // UTILITIES_H
