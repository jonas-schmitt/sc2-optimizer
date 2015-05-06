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
    int mShieldRegenCount = 0;
public:
//    ProtossUnit();

//    ProtossUnit(string name);

//    ProtossUnit(const UnitStats& baseStats);

//    ProtossUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

//    ProtossUnit(ProtossUnit const& protossUnit);

//    ProtossUnit(BaseUnit const& baseUnit);

    void subShield(double const value);

    void subHealth(double const value);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer <= 0)
        {
            if(mShieldRegenCount <= 0 && mStats.shield < mStats.maxShield && mStats.health > EPS)
            {
                BaseUnit::addShield(2.0 * state.regenerationUpdate / 1000);
            }
            else if(mShieldRegenCount > 0)
            {
                --mShieldRegenCount;
            }

        }
    }

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        regenerate(own);
    }

    void initUpgrades(vector<int> const& flags);
};




// Protoss Units

class Probe final : public ProtossUnit
{};

class Zealot final : public ProtossUnit
{
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
                if(enemy->attackPossible(*this))
                {
                    applyCharge = true;
                    break;
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
    void initUpgrades(vector<int> const& flags);
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(mChargeAvail)
        {
            charge(other.unitList);
        }
        ProtossUnit::timestep(own, other);
    }

};

class Stalker final : public ProtossUnit
{};

class Sentry final : public ProtossUnit
{};

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
