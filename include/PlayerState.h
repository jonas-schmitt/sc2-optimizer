#ifndef _PLAYERSTATE_H_
#define _PLAYERSTATE_H_

#include<list>
#include<string>
#include<utility>
#include<algorithm>
#include<functional>
#include<vector>

#include"Utilities.h"

using std::string;
using std::list;
using std::pair;
using std::function;
using std::vector;

template<typename Race> struct PotentialField : public Race
{
    PotentialField(Vec2D p, function<Vec2D(Vec2D const&,typename Race::BUT const&)> f)
        : pos(p), func(f)
    {}
    Vec2D computeForce(typename Race::BUT const& unit) const
    {
        return func(pos,unit);
    }

private:
    Vec2D pos;
    function<Vec2D(Vec2D const&,typename Race::BUT const&)> func;

};


struct SimulationResult
{
    double damage = 0;
    double damagePercent = 0;
    double minerals = 0, gas = 0;
    size_t survivors = 0;
    vector<vector<Vec2Df>> paths;
};

template <class Race> struct PlayerState : public Race
{

    Vec2D startPos;

    Vec2D minPos;
    Vec2D maxPos;
    Vec2D fieldSize;

    list<PotentialField<Race>> potentialList;

    typedef typename std::list<typename Race::RUT *>::iterator unitIterator;

    list<pair<typename Race::UT0, unitIterator>> unitList0;
    list<pair<typename Race::UT1, unitIterator>> unitList1;
    list<pair<typename Race::UT2, unitIterator>> unitList2;
    list<pair<typename Race::UT3, unitIterator>> unitList3;
    list<pair<typename Race::UT4, unitIterator>> unitList4;
    list<pair<typename Race::UT5, unitIterator>> unitList5;
    list<pair<typename Race::UT6, unitIterator>> unitList6;
    list<pair<typename Race::UT7, unitIterator>> unitList7;
    list<pair<typename Race::UT8, unitIterator>> unitList8;
    list<pair<typename Race::UT9, unitIterator>> unitList9;
    list<pair<typename Race::UT10, unitIterator>> unitList10;
    list<pair<typename Race::UT11, unitIterator>> unitList11;
    list<pair<typename Race::UT12, unitIterator>> unitList12;
    list<pair<typename Race::UT13, unitIterator>> unitList13;
    list<pair<typename Race::UT14, unitIterator>> unitList14;
    list<pair<typename Race::UT15, unitIterator>> unitList15;
    list<pair<typename Race::UT16, unitIterator>> unitList16;
    list<pair<typename Race::UT17, unitIterator>> unitList17;

    list<typename Race::RUT*>unitList;
    size_t unitCount;

    int timeSlice = 10;
    long time = 0;
    int regenerationTimer = 0;

    int regenerationUpdate = 1000;

    template<typename T> void timestep(PlayerState<T>& other)
    {
        for(auto& elem : this->unitList0)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList1)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList2)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList3)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList4)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList5)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList6)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList7)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList8)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList9)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList10)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList11)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList12)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList13)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList14)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList15)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList16)
        {
            elem.first.timestep(*this, other);
        }
        for(auto& elem : this->unitList17)
        {
            elem.first.timestep(*this, other);
        }
        if(regenerationTimer < 0)
        {
            regenerationTimer = regenerationUpdate;
        }
        time += timeSlice;
        regenerationTimer -= timeSlice;
    }


    template<typename T> void eraseUnit(typename list<pair<T,unitIterator>>::iterator it, list<pair<T,unitIterator>>& unitListT)
    {
        unitList.erase((*it).second);
        unitListT.erase(it);
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
            unit->resetHealth();
            unit->resetShield();
            unit->resetEnergy();
            unit->resetPos();
            unit->resetTimer();
        }
        time = 0;
        regenerationTimer = regenerationUpdate;

    }

    template<typename T>
    void removeZombies(list<pair<T,unitIterator>>& unitListT)
    {
        if(unitListT.empty())
        {
            return;
        }

        for(auto it = unitListT.begin(); it != unitListT.end(); ++it)
        {
            if(it->first.getHealth() < EPS)
            {
                auto delIt = it;
                ++it;
                eraseUnit(delIt,unitListT);
            }
        }
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
        auto cmp = [] (vector<Vec2Df> const& a, vector<Vec2Df> const& b)
        {
            return a.size() > b.size();
        };
        sort(res.paths.begin(),res.paths.end(),cmp);

        return res;
    }


};

#endif
