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
    array<double, 13> X;

    double mHash;

    void initHash()
    {
        mHash = 0.0;
        for(double const val : X)
        {
            mHash += val;
        }

    }

public:
    Fitness fitness;
    int count;

    // Default constructor, creates random genes
    UnitGenes()
    {
        std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> dist(MIN,MAX);
        for(int i = 0; i < 13; ++i)
        {
            X[i] = dist(gen);
        }

        count = 0;
        initHash();

    }
    UnitGenes(double x)
    {
        X.fill(x);
        count = 0;
        initHash();

    }

    UnitGenes(initializer_list<double> L)
    {
        if(L.size() != 13)
        {
            throw std::invalid_argument("UnitGenes::UnitGenes(initializer_list<double>): The initializer list must contain 12 arguments");
        }
        auto it = L.begin();
        for(int i = 0; i < 13; ++i)
        {
            if(*it > MAX)
            {
                throw std::invalid_argument("UnitGenes::UnitGenes(initializer_list<double>): The first 5 values must be smaller than 2^8-1");
            }
            X[i] = *it;
            ++it;
        }

        count = 0;
        initHash();
    }

    UnitGenes(vector<double> const& L)
    {
        if(L.size() != 13)
        {
            throw std::invalid_argument("UnitGenes::UnitGenes(vector<double>): The vector must contain 12 arguments");
        }
        auto it = L.begin();
        for(int i = 0; i < 13; ++i)
        {
            if(*it > MAX)
            {
                throw std::invalid_argument("UnitGenes::UnitGenes(vector<double>): The first 5 values must be smaller than 2^8-1");
            }
            X[i] = *it;
            ++it;
        }

        count = 0;
        initHash();
    }

    // copy constructor
    UnitGenes(UnitGenes const& clone)
    {
        for(int i = 0; i < 13; ++i)
        {
            X[i] = clone.get(i);
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
        for(int i = 0; i < 13; ++i)
        {
            int k = dist(gen);
            if(k == 0)
            {
                X[i] = mother.get(i);
            }
            else if(k == 1)
            {
                X[i] = father.get(i);
            }
            else
            {
                X[i] = mother.get(i)/2+father.get(i)/2;
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
        size_t number = static_cast<size_t>(std::min(mutationRate*13.,13.));
        unordered_set<size_t> positions;
        std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist1(0,12);
        while(positions.size() < number)
        {
            positions.insert(dist1(gen));
        }
        std::uniform_real_distribution<double> dist2(MIN,MAX);
        for(int i = 0; i < 13; ++i)
        {
            if(positions.count(i) == 1)
            {
                X[i] = dist2(gen);
            }
            else
            {
                X[i] = original.get(i);
            }
        }

        count = 0;
        initHash();
    }

    double get(double const pos) const
    {
        return X[pos];
    }

    void set(int const pos, double const x)
    {
        X[pos] = x;
    }

    double getHash() const
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
    double computeDistance(UnitGenes const& other) const
    {

        double dist = 0;
        for(int i = 0; i < 13; ++i)
        {
            dist += std::abs(this->get(i)-other.get(i));
        }
        return dist;
    }
};

inline std::ostream& operator<<(std::ostream& out, UnitGenes const& genes)
{
    out << "X:";
    for(int i = 0; i < 13; ++i)
    {
        out << '\t' << std::to_string(genes.get(i));
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
