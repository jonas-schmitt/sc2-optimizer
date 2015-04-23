#ifndef UTILITIES_H
#define UTILITIES_H
#include <cmath>

static double const LIMIT = 10000;
static double const EPS = 1e-11;
static double const GASTOMINERALS = 1.5;
static double const STDLEN = 1./std::sqrt(2);
//static long const INTLIMIT = std::numeric_limits<int>::max();
//static long const BOUND = std::pow(2,16);

static int const MIN = 0.0;
static int const MAX = 1.0; //static_cast<int>(std::min(INTLIMIT,BOUND));
static double const MAX_INV = 1.0/MAX;


struct Vec2D
{

    double x;
    double y;

    Vec2D() : x(0.0), y(0.0) {}
    Vec2D(double val) : x(val), y(val) {}
    Vec2D(double x_init, double y_init) : x(x_init), y(y_init) {}
    Vec2D(Vec2D const& vec) : x(vec.x), y(vec.y) {}

    double computeLength() const
    {
        return std::sqrt(x*x+y*y);
    }
    Vec2D getNormedVec(double const len) const
    {
        Vec2D result(x / len, y / len);
        if(std::isinf (result.x) || std::isnan(result.x) || std::isinf (result.y) || std::isnan(result.y) || len < EPS)
        {
            // choose direction that brings the unit as much away from the border as possible

            if(x > EPS)
            {
                result.x = STDLEN;
            }
            else if(x < -EPS)
            {
                result.x = -STDLEN;
            }
            else
            {
                result.x = 0;
            }
            if(y > EPS)
            {
                result.y = STDLEN;
            }
            else if(y < -EPS)
            {
                result.y = -STDLEN;
            }
            else
            {
                result.y = 0;
            }
        }
        return result;
    }
    Vec2D getNormedVec() const
    {
        return getNormedVec (computeLength ());
    }
};
struct Vec2Df
{

    float x;
    float y;

    Vec2Df() : x(0.0), y(0.0) {}
    Vec2Df(float val) : x(val), y(val) {}
    Vec2Df(float x_init, float y_init) : x(x_init), y(y_init) {}
    Vec2Df(Vec2Df const& vec) : x(vec.x), y(vec.y) {}
    Vec2Df(Vec2D const& vec) : x(static_cast<float>(vec.x)), y(static_cast<float>(vec.y)) {}

    float computeLength() const
    {
        return std::sqrt(x*x+y*y);
    }
    Vec2Df getNormedVec(float const len) const
    {
        Vec2Df result(x / len, y / len);
        if(std::isinf (result.x) || std::isnan(result.x) || std::isinf (result.y) || std::isnan(result.y) || len < EPS)
        {
            // choose direction that brings the unit as much away from the border as possible

            if(x > EPS)
            {
                result.x = STDLEN;
            }
            else if(x < -EPS)
            {
                result.x = -STDLEN;
            }
            else
            {
                result.x = 0;
            }
            if(y > EPS)
            {
                result.y = STDLEN;
            }
            else if(y < -EPS)
            {
                result.y = -STDLEN;
            }
            else
            {
                result.y = 0;
            }
        }
        return result;
    }
    Vec2Df getNormedVec() const
    {
        return getNormedVec(computeLength ());
    }
};


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
