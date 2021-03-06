
#include<utility>
#include<functional>
#include<cmath>
#include<list>
#include "../include/Terran.h"


void TerranBioUnit::initUpgrades (std::vector<int> const& flags)
{
    TerranUnit::initUpgrades(flags);
    if(flags.size() < 3)
    {
        return;
    }
    if(flags[2] == 1)
    {
        mStimpackAvail = true;
    }
}

void Marine::initUpgrades(std::vector<int> const &flags)
{
    TerranBioUnit::initUpgrades (flags);
    if(flags.size () < 4)
    {
        return;
    }
    if(flags[3] == 1)
    {
        mStats.maxHealth += 10.0;
        mStats.health += 10.0;
    }
}

void Marauder::concussiveShells ()
{
    if(mTarget != nullptr && !mTarget->mCSAffected)
    {
        mTarget->multSpeed (0.5);
        mTarget->mCSAffected = true;
        int const duration = 1500;
        mCSData.emplace_back(duration, mTarget);
        mTarget = nullptr;
    }
    if(!mCSData.empty())
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

void Marauder::initUpgrades(const std::vector<int> &flags)
{
    TerranBioUnit::initUpgrades (flags);
    if(flags.size () < 4)
    {
        return;
    }
    if(flags[3] == 1)
    {
        mCSAvail = true;
    }
}



//Reaper::Reaper()
//    : TerranUnit(), mHealthRegenCount(0)
//{}

//Reaper::Reaper(string name)
//    : TerranUnit(name), mHealthRegenCount(0)
//{}

//Reaper::Reaper(const UnitStats& baseStats)
//    : TerranUnit(baseStats), mHealthRegenCount(0)
//{}

//Reaper::Reaper(const UnitStats& baseStats, Vec2D min, Vec2D max)
//    : TerranUnit(baseStats,min,max), mHealthRegenCount(0)
//{}

//Reaper::Reaper(Reaper const& Reaper)
// : TerranUnit(Reaper), mHealthRegenCount(0)
//{}

//Reaper::Reaper(BaseUnit const& baseUnit)
//    : TerranUnit(baseUnit), mHealthRegenCount(0)
//{}

//Reaper::Reaper(TerranUnit const& terranUnit)
//    : TerranUnit(terranUnit), mHealthRegenCount(0)
//{}

void Reaper::subHealth(double const value)
{
    mHealthRegenCount = 10;
    TerranUnit::subHealth(value);
}

void Hellion::initUpgrades (std::vector<int> const& flags)
{
    TerranUnit::initUpgrades(flags);
    if(flags.size() < 3)
    {
        return;
    }
    if(flags[2] == 1)
    {
        for(Bonus& bonus : mStats.bonuses)
        {
            if(std::find(bonus.attributes.begin(), bonus.attributes.end(), Attribute::light) != bonus.attributes.end())
            {
                bonus.base += 5.0;
            }
        }
    }
}

void Hellbat::initUpgrades (std::vector<int> const& flags)
{
    TerranUnit::initUpgrades(flags);
    if(flags.size() < 3)
    {
        return;
    }
    if(flags[2] == 1)
    {
        for(Bonus& bonus : mStats.bonuses)
        {
            if(std::find(bonus.attributes.begin(), bonus.attributes.end(), Attribute::light) != bonus.attributes.end())
            {
                bonus.base += 12.0;
            }
        }
    }
}




