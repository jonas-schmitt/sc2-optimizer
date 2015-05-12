#include<utility>
#include<functional>
#include<cmath>
#include<list>

#include "../include/Protoss.h"





//ProtossUnit::ProtossUnit()
//    : BaseUnit(), mShieldRegenCount(0)
//{}

//ProtossUnit::ProtossUnit(string name)
//    : BaseUnit(name), mShieldRegenCount(0)
//{}

//ProtossUnit::ProtossUnit(const UnitStats& baseStats)
//    : BaseUnit(baseStats), mShieldRegenCount(0)
//{}

//ProtossUnit::ProtossUnit(const UnitStats& baseStats, Vec2D min, Vec2D max)
//    : BaseUnit(baseStats,min,max), mShieldRegenCount(0)
//{}

//ProtossUnit::ProtossUnit(ProtossUnit const& protossUnit)
// : BaseUnit(protossUnit), mShieldRegenCount(0)
//{}

//ProtossUnit::ProtossUnit(BaseUnit const& baseUnit)
//    : BaseUnit(baseUnit), mShieldRegenCount(0)
//{}


void ProtossUnit::subShield(double const value)
{
    mShieldRegenCount = 10;
    BaseUnit::subShield(value);
}

void ProtossUnit::subHealth(double const value)
{
    mShieldRegenCount = 10;
    BaseUnit::subHealth(value);
}

void ProtossUnit::initUpgrades (vector<int> const &flags)
{
    BaseUnit::initUpgrades(flags);
    if(flags.size() < 3)
    {
        return;
    }
    mShieldUpgrade = flags[2];
}

void Zealot::initUpgrades(const vector<int> &flags)
{
    ProtossUnit::initUpgrades(flags);
    if(flags.size() < 4)
    {
        return;
    }
    if(flags[3] == 1)
    {
        mChargeAvail = true;
        mStats.speed += 0.5;
    }

}

void Stalker::initUpgrades (const vector<int> &flags)
{
    ProtossUnit::initUpgrades(flags);
    if(flags.size() < 4)
    {
        return;
    }
    if(flags[3] == 1)
    {
        mBlinkAvail = true;
    }
}
