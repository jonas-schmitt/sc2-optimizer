#include<utility>
#include<functional>
#include<cmath>
#include<list>

#include "../include/Zerg.h"



//ZergUnit::ZergUnit()
//    : BaseUnit()
//{}

//ZergUnit::ZergUnit(string name)
//    : BaseUnit(name)
//{}

//ZergUnit::ZergUnit(const UnitStats& baseStats)
//    : BaseUnit(baseStats)
//{}

//ZergUnit::ZergUnit(const UnitStats& baseStats, Vec2D min, Vec2D max)
//    : BaseUnit(baseStats,min,max)
//{}

//ZergUnit::ZergUnit(ZergUnit const& zergUnit)
// : BaseUnit(zergUnit)
//{}

//ZergUnit::ZergUnit(BaseUnit const& baseUnit)
//    : BaseUnit(baseUnit)
//{}

void ZergUnit::initUpgrades (vector<int> const &flags)
{
    BaseUnit::initUpgrades (flags);
    if(flags.size () < 3)
    {
        return;
    }
    if(flags[2] == 1)
    {
        mStats.speed *= mStats.creepMultiplier;
    }
}

void Zergling::initUpgrades (vector<int> const& flags)
{
    ZergUnit::initUpgrades (flags);
    if(flags.size () < 5)
    {
        return;
    }
    if(flags[3] == 1)
    {
        mStats.speed *= 1.6;
        mMoveDist *= 1.6;
    }
    if(flags[4] == 1)
    {
        mStats.gaCooldown -= 109;
    }
}

void Baneling::initUpgrades (vector<int> const& flags)
{
    ZergUnit::initUpgrades (flags);
    if(flags.size () < 4)
    {
        return;
    }
    if(flags[3] == 1)
    {
        bool const creep = flags[2] == 1;
        if(creep) mStats.speed /= mStats.creepMultiplier;
        mStats.speed += 0.4531;
        if(creep) mStats.speed *= mStats.creepMultiplier;
    }


}



void Roach::initUpgrades (vector<int> const& flags)
{
    ZergUnit::initUpgrades (flags);
    if(flags.size () < 4)
    {
        return;
    }
    if(flags[3] == 1)
    {
        mStats.speed += 0.75;
    }

}


