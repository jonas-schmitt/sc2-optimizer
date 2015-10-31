#ifndef _PLAYERSTATE_H_
#define _PLAYERSTATE_H_


#include<utility>
#include<algorithm>
#include<functional>
#include<vector>
#include<queue>

#include"Utilities.h"
#include"Chromosome.h"


using std::vector;
using std::deque;

template<typename Race> struct PotentialField final : public Race
{
    PotentialField(Vec2D p, std::function<Vec2D(Vec2D const&,typename Race::BUT const&)> f)
        : pos(p), func(f)
    {}
    Vec2D computeForce(typename Race::BUT const& unit) const
    {
        return func(pos,unit);
    }

private:
    Vec2D pos;
    std::function<Vec2D(Vec2D const&,typename Race::BUT const&)> func;

};


struct SimulationResult final
{
    double damage = 0;
    double damagePercent = 0;
    double minerals = 0, gas = 0;
    size_t survivors = 0;
    std::vector<std::vector<Vec2Df>> paths;
};

template <class Race> struct PlayerState final : public Race
{

    bool forceFieldPlaced = false;
    Vec2D startPos;

    Vec2D minPos;
    Vec2D maxPos;
    Vec2D fieldSize;
    double maxUnitSize;

    std::vector<PotentialField<Race>> potentialList;
    std::deque<std::pair<int,PotentialField<Race>>> forceFieldQueue;

    Chromosome chromosome;


    std::vector<typename Race::UT0> unitList0;
    std::vector<typename Race::UT1> unitList1;
    std::vector<typename Race::UT2> unitList2;
    std::vector<typename Race::UT3> unitList3;
    std::vector<typename Race::UT4> unitList4;
    std::vector<typename Race::UT5> unitList5;
    std::vector<typename Race::UT6> unitList6;
    std::vector<typename Race::UT7> unitList7;
    std::vector<typename Race::UT8> unitList8;
    std::vector<typename Race::UT9> unitList9;
    std::vector<typename Race::UT10> unitList10;
    std::vector<typename Race::UT11> unitList11;
    std::vector<typename Race::UT12> unitList12;
    std::vector<typename Race::UT13> unitList13;
    std::vector<typename Race::UT14> unitList14;
    std::vector<typename Race::UT15> unitList15;
    std::vector<typename Race::UT16> unitList16;
    std::vector<typename Race::UT17> unitList17;

    std::vector<typename Race::RUT*>unitList;
    size_t unitCount;

    int timeSlice = 10;
    long time = 0;
    int regenerationTimer = 0;

    int regenerationUpdate = 1000;

    int NGenes = 0;




    template<typename T> void timestep(PlayerState<T>& other)
    {

        for(auto& elem : this->unitList0)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList1)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList2)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList3)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList4)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList5)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList6)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList7)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList8)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList9)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList10)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList11)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList12)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList13)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList14)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList15)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList16)
        {
            elem.timestep(*this, other);
        }
        for(auto& elem : this->unitList17)
        {
            elem.timestep(*this, other);
        }

        if(!forceFieldQueue.empty())
        {
            for(auto& forceField : forceFieldQueue)
            {
                forceField.first -= timeSlice;
            }
            if(forceFieldQueue.front().first <= 0)
            {
                forceFieldQueue.pop_front();
            }
        }
        if(regenerationTimer > 0)
        {
            regenerationTimer -= timeSlice;
        }
        else
        {
            regenerationTimer = regenerationUpdate;
        }
        time += timeSlice;

    }



    void clear()
    {
        unitList.clear();
        unitList0.clear();
        unitList1.clear();
        unitList2.clear();
        unitList3.clear();
        unitList4.clear();
        unitList5.clear();
        unitList6.clear();
        unitList7.clear();
        unitList8.clear();
        unitList9.clear();
        unitList10.clear();
        unitList11.clear();
        unitList12.clear();
        unitList13.clear();
        unitList14.clear();
        unitList15.clear();
        unitList16.clear();
        unitList17.clear();
        potentialList.clear();
    }
    void reset()
    {
        unitCount = unitList.size();
        for(auto unit : unitList)
        {
            unit->reset();
        }
        time = 0;
        regenerationTimer = regenerationUpdate;

    }



    SimulationResult calculateResult() const
    {
        SimulationResult res;
        double maxDamage = 0;
        res.paths.reserve(unitList.size());
        for(auto unit : unitList)
        {
            maxDamage += unit->getMaxHealth()+unit->getMaxShield();
            res.damage += unit->getMaxHealth()+unit->getMaxShield()-unit->getHealth()-unit->getShield();
            res.minerals += unit->getMinerals();
            res.gas += unit->getGas();
            if(unit->getHealth() > EPS)
            {
                ++res.survivors;
            }
            res.paths.push_back(unit->getPath());
        }
        res.damagePercent = res.damage*100. / maxDamage;
        auto cmp = [] (std::vector<Vec2Df> const& a, std::vector<Vec2Df> const& b)
        {
            return a.size() > b.size();
        };
        sort(res.paths.begin(),res.paths.end(),cmp);

        return res;
    }


};

#endif
