#ifndef TERRAN_H
#define TERRAN_H
#include<utility>
#include<functional>
#include<cmath>
#include<list>
#include<chrono>
#include<random>
#include<string>
#include<vector>
#include<typeinfo>
#include<array>
#include<deque>


#include "PlayerState.h"
#include "Utilities.h"
#include "UnitGenes.h"
#include "BaseUnit.h"

using std::pair;
using std::tuple;
using std::function;
using std::pow;
using std::sqrt;
using std::list;
using std::string;
using std::vector;
using std::array;
using std::deque;


class TerranUnit : public BaseUnit
{
public:
//    TerranUnit();

//    TerranUnit(string name);

//    TerranUnit(const UnitStats& baseStats);

//    TerranUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

//    TerranUnit(TerranUnit const& terranUnit);

//    TerranUnit(BaseUnit const& baseUnit);

    template <typename T> void regenerate(PlayerState<T>&) {}

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        regenerate(own);
    }
    //void initUpgrades(vector<int> const& flags);
};

class TerranBioUnit : public TerranUnit
{
    bool mStimpackAvail = false;
    bool mStimpack = false;
    int mStimpackTimer = 0;

    template<typename T> void stimpack(vector<T *> const & unitList)
    {
        if(mStimpackTimer <= 0)
        {
            if(mStimpack)
            {
                mStats.gaCooldown *= 2.0;
                mStats.aaCooldown *= 2.0;
                //mStats.speed *= 0.5;
                mMoveDist *= 0.5;
                mStimpack = false;
            }
            double const multiplier = 2.0;
            double const threshold = 20.0;

            bool applyStimpack = false;
            if(mStats.health > threshold)
            {
                for(auto const unit : unitList)
                {
                    if(computeDistance(*unit) < mMoveDist*multiplier + computeRange(*unit))
                    {
                        applyStimpack = true;
                        break;
                    }
                }
            }
            if(applyStimpack)
            {
                mStats.gaCooldown *= 0.5;
                mStats.aaCooldown *= 0.5;
                //mStats.speed *= 2.0;
                mMoveDist *= 2.0;
                mStats.health -= 10.0;
                mStimpackTimer = 15000;
                mStimpack = true;
                return;
            }
        }
        mStimpackTimer -= mTimeSlice;
    }

public:
    void initUpgrades(vector<int> const& flags);
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(mStimpackAvail)
        {
            stimpack(other.unitList);
        }
        TerranUnit::timestep(own, other);
    }
};


// classes for every specific Unit

// Terran Units

class SCV final : public TerranUnit
{};

class Marine final : public TerranBioUnit
{
public:
    void initUpgrades(vector<int> const& flags);
};

class Marauder final : public TerranBioUnit
{
private:
    bool mCSAvail = false;
    deque<pair<int,BaseUnit *>> mCSData;

    void concussiveShells();

public:
    void initUpgrades(vector<int> const& flags);
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        TerranBioUnit::timestep(own, other);
        if(mCSAvail)
        {
            concussiveShells();
            if(mCSData.size() > 0)
            {
                if(mCSData.front().first <= 0)
                {
                    BaseUnit& unit = *mCSData.front().second;
                    unit.multSpeed (2.0);
                    unit.mCSAffected = false;
                    mCSData.pop_front();
                }
                for(auto& el : mCSData)
                {
                    el.first -= mTimeSlice;
                }
            }
        }
    }
};

class Reaper final : public TerranUnit
{
protected:
    int mHealthRegenCount = 0;
public:

//    Reaper();

//    Reaper(string name);

//    Reaper(const UnitStats& baseStats);

//    Reaper(const UnitStats& baseStats, Vec2D min, Vec2D max);

//    Reaper(Reaper const& Reaper);

//    Reaper(BaseUnit const& baseUnit);

//    Reaper(TerranUnit const& terranUnit);

    void subHealth(double const value);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer <= 0)
        {
            if(mHealthRegenCount <= 0 && mStats.health < mStats.maxHealth && mStats.health > EPS)
            {
                TerranUnit::addHealth(2.0 * state.regenerationUpdate / 1000);
            }
            else
            {
                if(mHealthRegenCount > 0)
                {
                    --mHealthRegenCount;
                }
            }
        }
    }

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        TerranUnit::timestep(own, other);
        regenerate(own);
    }
};

class Ghost final : public TerranUnit
{};

class Hellion final : public TerranUnit
{
public:
    void initUpgrades(vector<int> const& flags);
};

class Hellbat final : public TerranUnit
{
public:
    void initUpgrades(vector<int> const& flags);
};

class SiegeTank final : public TerranUnit
{};

class WidowMine final : public TerranUnit
{};

class Thor final : public TerranUnit
{};

class Viking final : public TerranUnit
{};

class Medivac final : public TerranUnit
{};

class Raven final : public TerranUnit
{};

class Banshee final : public TerranUnit
{};

class Battlecruiser final : public TerranUnit
{};

#endif // TERRAN_H
