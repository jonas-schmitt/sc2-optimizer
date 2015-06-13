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
#include "Chromosome.h"
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
                addHealth(0.27 * 1e-3 * state.regenerationUpdate);
            }
        }
    }

public:

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
                addHealth(0.27 * 1e-3 * state.regenerationUpdate);
            }
            if(mStats.energy < mStats.maxEnergy)
            {
                mStats.energy += 0.5625 * 1e-3 * state.regenerationUpdate;
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
    int const NGENES = ZergUnit::NGENES + 1;
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        Queen::regenerate(own);
        Queen::transfuse(own.unitList);
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
    int const NGENES = ZergUnit::NGENES + 2;

    void initUpgrades (vector<int> const& flags);
    template<typename U, typename T> bool attack(PlayerState<U>& own, PlayerState<T>& other)
    {
        if(this->getHealth() < EPS || other.unitList.empty())
        {
            return false;
        }
        double sum = 0.0;
        double const threshold = this->getGroundAttack ();
        vector<pair<Damage,typename T::RUT *>> targets;
        for(auto enemy : other.unitList)
        {
            Damage damage;
            if(attackPossible (*enemy))
            {
                damage = computeDamage(*enemy);
                targets.emplace_back(damage, enemy);
                sum += damage.total;

            }
        }

        if(sum > threshold)
        {
            for(auto el : targets)
            {
                Damage const& damage = el.first;
                auto unit = el.second;
                if(damage.shield > 0)
                {
                    unit->subShield(damage.shield);
                }
                if(damage.health > 0)
                {
                    unit->subHealth(damage.health);
                }
                if(unit->getHealth() < EPS)
                {
                    --other.unitCount;
                }
            }
            this->setHealth(0.0);
            --own.unitCount;

            return true;
        }
        return false;
    }

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(!(Baneling::attack(own, other)))
        {
            ZergUnit::move(own, other);
        }
        ZergUnit::decMovementTimer();
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
