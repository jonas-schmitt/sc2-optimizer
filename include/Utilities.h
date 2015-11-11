#ifndef UTILITIES_H
#define UTILITIES_H
#include <cmath>
#include <cstddef>
#include <iostream>


static double constexpr LIMIT = 10000;
static double constexpr EPS = 1e-12;
static double constexpr GASTOMINERALS = 1.5;
static double constexpr STDLEN = 1./std::sqrt(2);
//static long const INTLIMIT = std::numeric_limits<int>::max();
//static long const BOUND = std::pow(2,16);

static double constexpr MIN = 0.0;
static double constexpr MAX = 1.0; //static_cast<int>(std::min(INTLIMIT,BOUND));
static double constexpr MAX_INV = 1.0/MAX;

static double constexpr INF = std::numeric_limits<double>::max();


struct Vec2D final
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

    double computeLengthSquared() const
    {
        return x*x+y*y;
    }

    Vec2D getNormedVec(double const len) const
    {
        /*if(len < EPS)
        {
            Vec2D result;
            if(x > 0.0)
            {
                result.x = STDLEN;
            }
            else if(x < 0.0)
            {
                result.x = -STDLEN;
            }
            if(y > 0.0)
            {
                result.y = STDLEN;
            }
            else if(y < 0.0)
            {
                result.y = -STDLEN;
            }
            return result;
        }*/
        double const tmp = 1.0/len;
        double const x_res = tmp * x;
        double const y_res = tmp * y;
        if(std::isinf (x_res) || std::isnan(x_res) || std::isinf (y_res) || std::isnan(y_res))
        {
            // choose direction that brings the unit as much away from the border as possible
            Vec2D result;
            if(x > 0.0)
            {
                result.x = STDLEN;
            }
            else if(x < 0.0)
            {
                result.x = -STDLEN;
            }
            if(y > 0.0)
            {
                result.y = STDLEN;
            }
            else if(y < 0.0)
            {
                result.y = -STDLEN;
            }
            return result;
        }
        return Vec2D(x_res, y_res);
    }
    Vec2D getNormedVec() const
    {
        return getNormedVec (computeLength ());
    }
};
struct Vec2Df final
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
        if(len < EPS)
        {
            Vec2Df result;
            if(x > 0.0)
            {
                result.x = STDLEN;
            }
            else if(x < 0.0)
            {
                result.x = -STDLEN;
            }
            if(y > 0.0)
            {
                result.y = STDLEN;
            }
            else if(y < 0.0)
            {
                result.y = -STDLEN;
            }
            return result;
        }
        float const tmp = 1.0/len;
        float const x_res = tmp * x;
        float const y_res = tmp * y;
        if(std::isinf (x_res) || std::isnan(x_res) || std::isinf (y_res) || std::isnan(y_res))
        {
            // choose direction that brings the unit as much away from the border as possible
            Vec2Df result;
            if(x > 0.0)
            {
                result.x = STDLEN;
            }
            else if(x < 0.0)
            {
                result.x = -STDLEN;
            }
            if(y > 0.0)
            {
                result.y = STDLEN;
            }
            else if(y < 0.0)
            {
                result.y = -STDLEN;
            }
            return result;
        }
        return Vec2Df(x_res, y_res);
    }
    Vec2Df getNormedVec() const
    {
        return getNormedVec(computeLength ());
    }
};

struct Statistics
{
    double mean;
    double max;
    double sum;
    double stdev;
    double onlinePerformance;
    double offlinePerformance;
    size_t iteration = 0;
    void print()
    {
        std::cout << "Total: " << sum << std::endl;
        std::cout << "Average: " << mean*100.0 << " %" << std::endl;
        std::cout << "Maximum: " << max*100.0 << " %" << std::endl;
        std::cout << "Standard Deviation: " << stdev*100 << " %" << std::endl;
        std::cout << "Online Performance: " << onlinePerformance*100.0/iteration << " %" << std::endl;
        std::cout << "Offline Performance: " << offlinePerformance*100.0/iteration << " %" << std::endl;
    }
};

enum class Player
{
    first, second
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


#endif // UTILITIES_H
