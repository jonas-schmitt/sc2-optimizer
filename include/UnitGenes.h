#ifndef UNITGENES_H
#define UNITGENES_H
#include<array>
#include<vector>
#include<unordered_set>
#include<random>
#include<chrono>
#include <iostream>
#include <initializer_list>

using std::array;
using std::unordered_set;
using std::vector;
using std::initializer_list;

struct UnitGenes
{
private:
    array<int, 5> X;
    array<int, 7> Y;
    size_t mHash;

    void initHash()
    {
        mHash = 0;
        for(int const value : X)
        {
            mHash ^= value + 0x9e3779b9 + (mHash<<6) + (mHash>>2);
        }
        for(int const value : Y)
        {
            mHash ^= value + 0x9e3779b9 + (mHash<<6) + (mHash>>2);
        }
    }

public:
    Fitness fitness;
    size_t count;

    // Default constructor, creates random genes
    UnitGenes()
    {
        std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist1(XMIN,XMAX);
        for(size_t i = 0; i < 5; ++i)
        {
            X[i] = dist1(gen);
        }

        std::uniform_int_distribution<int> dist2(YMIN,YMAX);
        for(size_t i = 0; i < 7; ++i)
        {
            Y[i] = dist2(gen);
        }
        count = 0;
        initHash();

    }
    UnitGenes(initializer_list<int> L)
    {
        if(L.size() != 12)
        {
            throw std::invalid_argument("UnitGenes::UnitGenes(initializer_list<int>): The initializer list must contain 12 arguments");
        }
        auto it = L.begin();
        for(size_t i = 0; i < 5; ++i)
        {
            if(*it > XMAX)
            {
                throw std::invalid_argument("UnitGenes::UnitGenes(initializer_list<int>): The first 5 values must be smaller than 2^8-1");
            }
            X[i] = *it;
            ++it;
        }
        for(size_t i = 0; i < 7; ++i)
        {
            if(*it > YMAX)
            {
                throw std::invalid_argument("UnitGenes::UnitGenes(initializer_list<int>): The last 7 values must be smaller than 2^16-1");
            }
            Y[i] = *it;
            ++it;
        }
        count = 0;
        initHash();
    }

    UnitGenes(vector<int> const& L)
    {
        if(L.size() != 12)
        {
            throw std::invalid_argument("UnitGenes::UnitGenes(vector<int>): The vector must contain 12 arguments");
        }
        auto it = L.begin();
        for(size_t i = 0; i < 5; ++i)
        {
            if(*it > XMAX)
            {
                throw std::invalid_argument("UnitGenes::UnitGenes(vector<int>): The first 5 values must be smaller than 2^8-1");
            }
            X[i] = *it;
            ++it;
        }
        for(size_t i = 0; i < 7; ++i)
        {
            if(*it > YMAX)
            {
                throw std::invalid_argument("UnitGenes::UnitGenes(vector<int>): The last 7 values must be smaller than 2^16-1");
            }
            Y[i] = *it;
            ++it;
        }
        count = 0;
        initHash();
    }

    // copy constructor
    UnitGenes(UnitGenes const& clone)
    {
        for(size_t i = 0; i < 5; ++i)
        {
            X[i] = clone.getX(i);
        }

        for(size_t i = 0; i < 7; ++i)
        {
            Y[i] = clone.getY(i);
        }
        fitness = clone.fitness;
        count = clone.count;
        mHash = clone.getHash();
    }
    // constructor implementing crossover
    UnitGenes(UnitGenes const& mother, UnitGenes const& father)
    {
        std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist(0,2);
        for(int i = 0; i < 5; ++i)
        {
            int k = dist(gen);
            if(k == 0)
            {
                X[i] = mother.getX(i);
            }
            else if(k == 1)
            {
                X[i] = father.getX(i);
            }
            else
            {
                X[i] = mother.getX(i)/2+father.getX(i)/2;
            }
        }
        for(size_t i = 0; i < 7; ++i)
        {
            int k = dist(gen);
            if(k == 0)
            {
                Y[i] = mother.getY(i);
            }
            else if(k == 1)
            {
                Y[i] = father.getY(i);
            }
            else
            {
                Y[i] = mother.getY(i)/2+father.getY(i)/2;
            }
        }
        count = 0;
        initHash();
    }


    // constructor implementing mutation
    UnitGenes(UnitGenes const& original, float const mutationRate)
    {
        if(mutationRate < 0+EPS || mutationRate > 1-EPS)
        {
            throw std::invalid_argument("UnitGenes::UnitGenes(UnitGenes const & mutant, float rate): The rate must be a value in (0,1)");
        }
        size_t number = static_cast<size_t>(std::min(mutationRate*12.,12.));
        unordered_set<size_t> positions;
        std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<size_t> dist1(0,11);
        while(positions.size() < number)
        {
            positions.insert(dist1(gen));
        }
        std::uniform_int_distribution<int> dist2(XMIN,XMAX);
        for(int i = 0; i < 5; ++i)
        {
            if(positions.count(i) == 1)
            {
                X[i] = dist2(gen);
            }
            else
            {
                X[i] = original.getX(i);
            }
        }
        std::uniform_int_distribution<int> dist3(YMIN,YMAX);
        for(int i = 0; i < 7; ++i)
        {
            if(positions.count(i+5) == 1)
            {
                Y[i] = dist3(gen);
            }
            else
            {
                Y[i] = original.getY(i);
            }
        }
        count = 0;
        initHash();
    }

    int getX(size_t const pos) const
    {
        return X[pos];
    }
    int getY(size_t const pos) const
    {
        return Y[pos];
    }
    void setX(size_t const pos, int const x)
    {
        X[pos] = x;
    }
    void setY(size_t const pos, int const y)
    {
        Y[pos] = y;
    }
    size_t getHash() const
    {
        return mHash;
    }

    bool operator<(UnitGenes const& genes) const
    {
        if(this->fitness.score != genes.fitness.score)
        {
            return this->fitness.score > genes.fitness.score;
        }
        else
        {
            return this->count > genes.count;
        }
    }
    long computeDistance(UnitGenes const& other) const
    {

        long distX = 0, distY = 0;
        for(size_t i = 0; i < 5; ++i)
        {
            distX += std::abs(this->getX(i)-other.getX(i));
            distY += std::abs(this->getY(i)-other.getY(i));
        }
        distY += std::abs(this->getY(5)-other.getY(5)) + std::abs(this->getY(6)-other.getY(6));
        return distX + distY;
    }
};

inline std::ostream& operator<<(std::ostream& out, UnitGenes const& genes)
{
    out << "X:";
    for(size_t i = 0; i < 5; ++i)
    {
        out << '\t' << std::to_string(genes.getX(i));
    }
    out << std::endl;
    out << "Y:";
    for(size_t i = 0; i < 7; ++i)
    {
        out << '\t' << std::to_string(genes.getY(i));
    }
    out << std::endl;
    out << "Score: " << genes.fitness.score << std::endl;
    out << "Damage: " << genes.fitness.damage*100. << '%' << std::endl;
    out << "Health: " << genes.fitness.health*100. << '%' << std::endl;
    out << "Count: " << genes.count << std::endl;
    return out;
}
inline std::ostream& operator<<(std::ostream& out, vector<UnitGenes> const& population)
{
    for (auto genes : population)
        out << genes << std::endl;
    return out;
}
#endif // UNITGENES_H
