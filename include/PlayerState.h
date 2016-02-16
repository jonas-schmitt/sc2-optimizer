#ifndef _PLAYERSTATE_H_
#define _PLAYERSTATE_H_


#include<utility>
#include<algorithm>
#include<functional>
#include<vector>
#include<queue>

#include"Utilities.h"
#include"Chromosome.h"


// Struct representing a potential field
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



// Struct containing all the data associated with a player
template <class Race> struct PlayerState final : public Race
{

    bool forceFieldPlaced = false;
    Vec2D startPos;

    Vec2D minPos;
    Vec2D maxPos;
    Vec2D fieldSize;
    double maxUnitSize;

    // Potential fields that are present independent of the units
    std::vector<PotentialField<Race>> potentialList;

    // Force fields (special ability of the sentry)
    // In case that no sentry is present, this data structure is always empty
    std::deque<std::pair<int,PotentialField<Race>>> forceFieldQueue;

    // Combined (concatenated) chromosome of all units
    Chromosome chromosome;


    // Vectors used for storing the actual units
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

    // Pointers to the base class of the race
    // Allows easy iteration over all units when no unit specific abilities and characteristics are required
    std::vector<typename Race::RUT*>unitList;

    // Number of units that are alive
    size_t unitCount;

    int timeSlice = 10;
    long time = 0;
    int regenerationTimer = 0;

    int regenerationUpdate = 1000;

    // Length of the combined chromosome
    int NGenes = 0;

    double APM = 400;

    bool kill = false;

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
        if(kill)
        {
            adjustActionsPerUnit();
            kill = false;
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
        kill = false;
        for(auto unit : unitList)
        {
            unit->reset();
        }
        time = 0;
        regenerationTimer = regenerationUpdate;

    }

    void decUnitCount()
    {
        --unitCount;
        kill = true;
    }

    void adjustActionsPerUnit()
    {
        double const actionsPerUnit = APM/unitCount;
        int const update = static_cast<int>(std::round(60000.0/actionsPerUnit));
        for(auto unit : unitList)
        {
            unit->setMovementUpdate(update);
        }
    }

};

#endif
