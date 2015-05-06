#ifndef ZERG_H
#define ZERG_H

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

class ZergUnit : public BaseUnit
{
protected:
    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer <= 0)
        {
            if(mStats.health > EPS && mStats.health < mStats.maxHealth)
            {
                addHealth(0.27 * state.regenerationUpdate * 1e-3);
            }
        }
    }

public:
//    ZergUnit();

//    ZergUnit(string name);

//    ZergUnit(const UnitStats& baseStats);

//    ZergUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

//    ZergUnit(ZergUnit const& zergUnit);

//    ZergUnit(BaseUnit const& baseUnit);

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        regenerate(own);
    }
    void initUpgrades(vector<int> const& flags);

};

// Zerg Units

class Drone final : public ZergUnit
{};

class Queen final : public ZergUnit
{
private:

    int mTransfusionTimer = 0;

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer <= 0)
        {
            if(mStats.health > EPS && mStats.health < mStats.maxHealth)
            {
                addHealth(0.27 * state.regenerationUpdate * 1e-3);
            }
            if(mStats.energy < mStats.maxEnergy)
            {
                mStats.energy += 0.5625 *state.regenerationUpdate * 1e-3;
                if(mStats.energy > mStats.maxEnergy)
                {
                    mStats.energy = mStats.maxEnergy;
                }
            }
        }
    }
    template<typename T> void transfuse(vector<T *>& unitList)
    {
        if(mTransfusionTimer <= 0 && mStats.energy >= 50)
        {
            T *target = nullptr;
            double max = 0.0;
            double const threshold = 100.0;
            for(T* unit : unitList)
            {
                if(unit == static_cast<T *>(this) || unit->getHealth() < EPS)
                {
                    continue;
                }
                double const lostHealth = unit->getMaxHealth() - unit->getHealth();
                if(lostHealth > threshold && computeDistance(*unit) - 7.0 < EPS)
                {
                    if(lostHealth > max)
                    {
                        max = lostHealth;
                        target = unit;
                    }
                }
            }
            if(max > 0.0)
            {
                target->incHealth(125.0);
                mStats.energy -= 50.0;
                mTransfusionTimer = 1000;
                return;
            }

        }
        if(mTransfusionTimer > 0)
        {
            mTransfusionTimer -= mTimeSlice;
        }

    }

public:

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        regenerate(own);
        transfuse(own.unitList);
    }
};

class Zergling final : public ZergUnit
{
public:
    void initUpgrades(vector<int> const& flags);
};

class Baneling final : public ZergUnit
{
public:
    void initUpgrades (vector<int> const& flags);
    template<typename T> bool attack(PlayerState<T>& other)
    {
        if(this->getHealth() < EPS || other.unitList.empty())
        {
            return false;
        }
        bool hasTarget = false;
        double const threshold = 0.5*this->getGroundAttack ();
        vector<pair<Damage,typename T::RUT *>> targets;
        for(auto enemy : other.unitList)
        {
            Damage damage;
            if(attackPossible (*enemy))
            {
                damage = computeDamage(*enemy);
                targets.emplace_back(damage, enemy);
                if(damage.total > threshold)
                {
                    hasTarget = true;
                }
            }
        }
        if(hasTarget)
        {
            for(auto el : targets)
            {
                Damage const& damage = el.first;
                auto unit = el.second;
                if(damage.shield > 0)
                {
                    unit->subShield(damage.shield);
                }
                unit->subHealth(damage.health);
            }
            this->setHealth(0.0);
            return true;
        }
        return false;
    }

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(!(attack(other)))
        {
            move(own, other);
        }
        decMovementTimer();
    }
};

class Roach final : public ZergUnit
{
public:
    void initUpgrades (vector<int> const& flags);
};

class Hydralisk final : public ZergUnit
{};

class Infestor final : public ZergUnit
{};

class SwarmHost final : public ZergUnit
{};

class Ultralisk :public ZergUnit
{};

class Overseer final : public ZergUnit
{};

class Mutalisk final : public ZergUnit
{};

class Corruptor final : public ZergUnit
{};

class BroodLord final : public ZergUnit
{};

class Viper final : public ZergUnit
{};

#endif // ZERG_H
