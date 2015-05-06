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
//public:
//    TerranUnit();

//    TerranUnit(string name);

//    TerranUnit(const UnitStats& baseStats);

//    TerranUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

//    TerranUnit(TerranUnit const& terranUnit);

//    TerranUnit(BaseUnit const& baseUnit);

//    template <typename T> void regenerate(PlayerState<T>&) {}

//    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
//    {
//        BaseUnit::timestep(own, other);
//        regenerate(own);
//    }
    //void initUpgrades(vector<int> const& flags);
};

class TerranBioUnit : public TerranUnit
{
private:
    bool mStimpackAvail = false;
    bool mStimpackActive = false;
    int mStimpackDurationTimer = 0;

    template<typename T> void stimpack(vector<T *> const & unitList)
    {
        if(mStimpackDurationTimer <= 0)
        {
            if(mStimpackActive)
            {
                mStats.gaCooldown *= 2.0;
                mStats.aaCooldown *= 2.0;
                //mStats.speed *= 0.5;
                mMoveDist *= 0.5;
                mStimpackActive = false;
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
                mStimpackDurationTimer = 15000;
                mStimpackActive = true;
                return;
            }
        }
        mStimpackDurationTimer -= mTimeSlice;
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
private:
    int mHealthRegenCount = 0;
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
public:

//    Reaper();

//    Reaper(string name);

//    Reaper(const UnitStats& baseStats);

//    Reaper(const UnitStats& baseStats, Vec2D min, Vec2D max);

//    Reaper(Reaper const& Reaper);

//    Reaper(BaseUnit const& baseUnit);

//    Reaper(TerranUnit const& terranUnit);

    void subHealth(double const value);

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
    template<typename T> bool attack(PlayerState<T>& other)
    {
        if(this->mAttackTimer > 0 || this->getHealth() < EPS || other.unitList.empty())
        {
            return false;
        }

        typename T::RUT *mainTarget = nullptr;
        double maxDamage = 0.0;
        bool kill = false;
        for(auto enemy : other.unitList)
        {
            Damage damage;
            if(attackPossible (*enemy))
            {
                damage = computeDamage(*enemy);
            }
            else
            {
                continue;
            }

            bool const newKill = 0.75 * damage.total > enemy->getHealth() + enemy->getShield();
            bool const higher = damage.total > maxDamage;
            if(newKill)
            {
                if(kill)
                {
                    if(higher)
                    {
                        maxDamage = damage.total;
                        mainTarget = enemy;
                        mPossibleDamage[enemy->getIdentifier()] = damage;
                    }
                }
                else
                {
                    maxDamage = damage.total;
                    mainTarget = enemy;
                    mPossibleDamage[enemy->getIdentifier()] = damage;
                }
            }
            else
            {
                if(higher)
                {
                    maxDamage = damage.total;
                    mainTarget = enemy;
                    mPossibleDamage[enemy->getIdentifier()] = damage;
                }
            }
        }
        if(maxDamage < EPS)
        {
            return false;
        }
        Vec2D distVec = computeDistance(*mainTarget);
        distVec = distVec.getNormedVec ();
        double const x1 = mPos.x;
        double const y1 = mPos.y;
        double const x2 = 5.0 * distVec.x + x1;
        double const y2 = 5.0 * distVec.y + y1;
        double const a = x2 - x1;
        double const b = y2 - y1;
        double const c = x2*y1 - y2*x1;
        double const d = 1.0/std::sqrt(a*a + b*b);

        vector<typename T::RUT *> targets;
        for(typename T::RUT *enemy : other.unitList)
        {
            if(mainTarget == enemy)
            {
                continue;
            }
            double const x0 = enemy->getX();
            double const y0 = enemy->getY();
            double const distFromLine = std::abs(b*x0 - a*y0 + c)*d;
            if(distFromLine - 0.15 < EPS)
            {
                double const delta_x1 = x0 - x1;
                double const delta_y1 = y0 - y1;
                double const delta_x2 = x0 - x2;
                double const delta_y2 = y0 - y2;
                double const dist1 = std::sqrt(delta_x1*delta_x1 + delta_y1*delta_y1);
                double const dist2 = std::sqrt(delta_x2*delta_x2 + delta_y2*delta_y2);
                if(dist1 - 5.15 < EPS && dist2 - 5.15 < EPS)
                {
                    targets.push_back(enemy);
                }
            }

        }

        if(BaseUnit::attack(*mainTarget))
        {
            --other.unitCount;
        }
        for(typename T::RUT *target : targets)
        {
            if(BaseUnit::attack(*target))
            {
                --other.unitCount;
            }
        }
        return true;
    }

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(!(attack(other)))
        {
            move(own, other);
        }
        decMovementTimer();
        decAttackTimer ();
    }
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
