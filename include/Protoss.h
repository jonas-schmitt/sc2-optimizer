#ifndef PROTOSS_H
#define PROTOSS_H
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

class ProtossUnit : public BaseUnit
{
protected:
    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer <= 0)
        {
            if(mShieldRegenCount <= 0 && mStats.shield < mStats.maxShield && mStats.health > EPS)
            {
                BaseUnit::addShield(2.0 * 1e-3 * state.regenerationUpdate);
            }
            else if(mShieldRegenCount > 0)
            {
                --mShieldRegenCount;
            }

        }
    }
    int mShieldRegenCount = 0;
public:

    void subShield(double const value);

    void subHealth(double const value);


    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        ProtossUnit::regenerate(own);
    }

    void initUpgrades(vector<int> const& flags);
};




// Protoss Units

class Probe final : public ProtossUnit
{};

class Zealot final : public ProtossUnit
{
private:
    bool mChargeAvail = false;
    bool mChargeActive = false;
    int mChargeAvailTimer = 0;
    int mChargeDurationTimer = 0;
    double const mApplyChargeFactor = 2.2;
    double const mRemoveChargeFactor = 1.0/mApplyChargeFactor;

    template<typename T> void charge(vector<T *> const & unitList)
    {

        if(mChargeAvailTimer <= 0 && !mChargeActive)
        {
            bool applyCharge = false;
            for(T *enemy : unitList)
            {
                if(!enemy->isDead())
                {

                    double const x = this->getX() - enemy->getX();
                    double const y = this->getY() - enemy->getY();
                    double dist = std::sqrt(x*x+y*y);
                    if(std::isnan(dist))
                    {
                        dist = 0.0;
                    }
                    double const threshold = enemy->computeRange(*this);
                    if(dist < threshold)
                    {
                        applyCharge = true;
                        break;
                    }
                }
            }
            if(applyCharge)
            {
                mStats.speed *= mApplyChargeFactor;
                mMoveDist *= mApplyChargeFactor;

                mChargeActive = true;
                mChargeDurationTimer = 3500;

                mChargeAvailTimer = 10000;
            }
        }
        if(mChargeAvailTimer > 0)
        {
            mChargeAvailTimer -= mTimeSlice;
        }

        if(mChargeDurationTimer > 0)
        {
            mChargeDurationTimer -= mTimeSlice;
        }
        else if(mChargeActive)
        {
            mChargeActive = false;
            mStats.speed *= mRemoveChargeFactor;
            mMoveDist *= mRemoveChargeFactor;
        }
    }

public:
    int const NGENES = ProtossUnit::NGENES + 1;
    void initUpgrades(vector<int> const& flags);
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(mChargeAvail)
        {
            Zealot::charge(other.unitList);
        }
        ProtossUnit::timestep(own, other);
    }

};

class Stalker final : public ProtossUnit
{
private:
    int mBlinkTimer = 0;
    bool mBlinkAvail = false;

    template<typename T> void blink(vector<T *> const & unitList)
    {
        int const threshold = 0;

        if(mBlinkTimer <= 0 && mStats.shield <= threshold)
        {
            double center_x = 0.0;
            double center_y = 0.0;
            for(T *enemy : unitList)
            {
                if(enemy->isDead())
                {
                    continue;
                }
                center_x += enemy->getX();
                center_y += enemy->getY();
            }
            double const n = 1.0/static_cast<double>(unitList.size());
            center_x *= n;
            center_y *= n;
            Vec2D distVec(center_x - mPos.x, center_y - mPos.y);
            distVec = std::move(distVec.getNormedVec());
            setPos(-8.0 * distVec.x + mPos.x, -8.0 * distVec.y + mPos.y);
            mBlinkTimer = 10000;
        }
        if(mBlinkTimer > 0)
        {
            mBlinkTimer -= mTimeSlice;
        }
    }

public:
    int const NGENES = ProtossUnit::NGENES + 1;
    void initUpgrades(vector<int> const& flags);
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(mBlinkAvail)
        {
            Stalker::blink(other.unitList);
        }
        ProtossUnit::timestep(own, other);
    }
};

class Sentry final : public ProtossUnit
{
private:
    int mForceFieldTimer = 0;
    bool mForceFieldPlaced = false;

    std::function<Vec2D(Vec2D const& pos, BaseUnit const& unit)> forceFieldFunc = [](Vec2D const& pos, BaseUnit const& unit)
    {
        double const radius = 1.7;
        Vec2D distVec(unit.getX() - pos.x, unit.getY() - pos.y);

        double const dist = distVec.computeLength();
        if(dist - radius < EPS)
        {
            distVec = std::move(distVec.getNormedVec(dist));
            return Vec2D(-1e5 * distVec.x, -1e5 * distVec.y);
        }
        return Vec2D(0.0);
    };

    template <typename T, typename U> void forceField(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(mForceFieldPlaced)
        {
            mForceFieldPlaced = false;
            other.forceFieldPlaced = false;
        }
        if(!other.forceFieldPlaced && mForceFieldTimer <= 0 && mStats.energy >= 50.0)
        {
            int enemyMovementTimer = mMovementUpdate;
            for(typename U::RUT *enemy : other.unitList)
            {
                if(enemy->getHealth() > EPS)
                {
                    enemyMovementTimer = enemy->getMovementTimer();
                    break;
                }
            }
            if(enemyMovementTimer > 0 && enemyMovementTimer <= mTimeSlice)
            {
                // apply force field
                Vec2D center(0.0);
                for(typename U::RUT *enemy : other.unitList)
                {
                    if(enemy->isDead())
                    {
                        continue;
                    }
                    center.x += enemy->getX();
                    center.y += enemy->getY();
                }
                double const n = 1.0/other.unitCount;
                center.x *= n;
                center.y *= n;
                Vec2D distVec(center.x - mPos.x, center.y - mPos.y);
                double const dist = distVec.computeLength();
                double const forceFieldRange = 9.0;
                double const forceFieldRadius = 1.7;
                double const threshold = 2.0 * (forceFieldRange + forceFieldRadius);
                if(dist < threshold)
                {
                    if(dist > forceFieldRange)
                    {
                        distVec = std::move(distVec.getNormedVec(dist));
                        center.x = mPos.x + forceFieldRange * distVec.x;
                        center.y = mPos.y + forceFieldRange * distVec.y;
                    }
                    // create Force Field

                    mForceFieldTimer = 15000;
                    own.forceFieldQueue.emplace_back(mForceFieldTimer, PotentialField<T>(center,forceFieldFunc));
                    other.forceFieldQueue.emplace_back(mForceFieldTimer, PotentialField<U>(center,forceFieldFunc));
                    mStats.energy -= 50.0;
                    mForceFieldPlaced = true;
                    other.forceFieldPlaced = true;
                }
            }
        }
        if(mForceFieldTimer > 0)
        {
            mForceFieldTimer -= mTimeSlice;
        }

    }

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer <= 0)
        {
            if(mShieldRegenCount <= 0 && mStats.shield < mStats.maxShield && mStats.health > EPS)
            {
                BaseUnit::addShield(2.0 * 1e-3 * state.regenerationUpdate);
            }
            else if(mShieldRegenCount > 0)
            {
                --mShieldRegenCount;
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

public:
    int const NGENES = ProtossUnit::NGENES + 1;
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        Sentry::forceField(own, other);
        BaseUnit::timestep(own, other);
        Sentry::regenerate(own);
    }
};

class HighTemplar final : public ProtossUnit
{};

class DarkTemplar final : public ProtossUnit
{};

class Immortal final : public ProtossUnit
{};

class Colossus final : public ProtossUnit
{};

class Archon final : public ProtossUnit
{};

class Observer final : public ProtossUnit
{};

class WarpPrism final : public ProtossUnit
{};

class Phoenix final : public ProtossUnit
{};

class VoidRay final : public ProtossUnit
{};

class Oracle final : public ProtossUnit
{};

class Carrier final : public ProtossUnit
{};

class Tempest final : public ProtossUnit
{};

class MothershipCore final : public ProtossUnit
{};

class Mothership final : public ProtossUnit
{};


#endif // PROTOSS_H
