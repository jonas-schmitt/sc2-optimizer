#include "../include/MicroSimulation.h"

#include<string>
#include<vector>
#include<list>
#include<memory>
#include<utility>
#include<algorithm>
#include<future>
#include<cmath>
#include<random>

using std::string;
using std::vector;
using std::shared_ptr;
using std::list;
using std::pair;

template <class T, class U>
MicroSimulation<T, U>::MicroSimulation(const Vec2D minPos, const Vec2D maxPos, const string& filePath1, const string& filePath2)
    : mMinPos(minPos), mMaxPos(maxPos), mFilePath1(filePath1), mFilePath2(filePath2), init1(filePath1), init2(filePath2)
{
    pl1.minPos = minPos;
    pl1.maxPos = maxPos;
    pl2.minPos = minPos;
    pl2.maxPos = maxPos;
}

template <class T, class U>
MicroSimulation<T,U>::MicroSimulation(MicroSimulation const &microSim)
    : MicroSimulation(microSim.getMinPos(), microSim.getMaxPos(), microSim.getFilePath1(), microSim.getFilePath2())
{}

template <class T, class U>
Vec2D MicroSimulation<T, U>::getMinPos() const
{
    return mMinPos;
}
template <class T, class U>
Vec2D MicroSimulation<T, U>::getMaxPos() const
{
    return mMaxPos;
}

template <class T, class U>
string MicroSimulation<T, U>::getFilePath1() const
{
    return mFilePath1;
}

template <class T, class U>
string MicroSimulation<T, U>::getFilePath2() const
{
    return mFilePath2;
}

template <class T, class U>
void MicroSimulation<T, U>::initPlayer1(const vector<string>& unitList)
{
    init1.init(unitList, mFilePath1, pl1);
    double const fieldSizeX = pl1.maxPos.x-pl1.minPos.x;
    double const fieldSizeY = pl1.maxPos.y-pl1.minPos.y;
    Vec2D startPos = Vec2D(pl1.minPos.x + fieldSizeX/20., pl1.minPos.y);
    for(auto unit : pl1.unitList)
    {
        unit->setMinPos(pl1.minPos);
        unit->setMaxPos(pl1.maxPos);
        startPos.y += fieldSizeY/pl1.unitList.size();
        unit->setStartPos(startPos);
        unit->resetPos();
    }
}

template <class T, class U>
void MicroSimulation<T, U>::initPlayer2(const vector<string>& unitList)
{
    init2.init(unitList, mFilePath2, pl2);
    double const fieldSizeX = pl2.maxPos.x-pl2.minPos.x;
    double const fieldSizeY = pl2.maxPos.y-pl2.minPos.y;
    Vec2D startPos = Vec2D(pl2.maxPos.x - fieldSizeX/20., pl2.minPos.y);
    for(auto unit : pl2.unitList)
    {
        unit->setMinPos(pl2.minPos);
        unit->setMaxPos(pl2.maxPos);
        startPos.y += fieldSizeY/pl2.unitList.size();
        unit->setStartPos(startPos);
        unit->resetPos();
    }
}



template <class T, class U>
void MicroSimulation<T, U>::initBothPlayers(const vector<string>& unitList1, const vector<string>& unitList2)
{
    initPlayer1 (unitList1);
    initPlayer2 (unitList2);

}

template<class T, class U>
template<typename V>
void MicroSimulation<T,U>::setPlayerChromosome(PlayerState<V>& pl, Chromosome const& chromosome)
{
    pl.chromosome = chromosome;
    size_t pos = 0;
    for(auto unitPtr : pl.unitList)
    {
        unitPtr->setChromosome(pl.chromosome);
    }
}

template<class T, class U>
void MicroSimulation<T,U>::setPlayer1Chromosome(Chromosome const& chromosome)
{
    setPlayerChromosome(pl1, chromosome);
}

template<class T, class U>
void MicroSimulation<T,U>::setPlayer2Chromosome(Chromosome const& chromosome)
{
    setPlayerChromosome(pl2, chromosome);
}

template<class T, class U>
void MicroSimulation<T,U>::clearPlayer1()
{
    pl1.clear();
}
template<class T, class U>
void MicroSimulation<T,U>::clearPlayer2()
{
    pl2.clear();
}


template<class T, class U>
void MicroSimulation<T,U>::clearBothPlayers()
{
    pl1.clear();
    pl2.clear();
}

template<class T, class U>
void MicroSimulation<T,U>::resetPlayer1()
{
    pl1.reset();
}

template<class T, class U>
void MicroSimulation<T,U>::resetPlayer2()
{
    pl2.reset();
}

template<class T, class U>
void MicroSimulation<T,U>::resetBothPlayers()
{
    pl1.reset();
    pl2.reset();
}



template <class T, class U>
void MicroSimulation<T, U>::setPlayer1Pos(Vec2D const pos)
{
    pl1.startPos = pos;
}

template <class T, class U>
void MicroSimulation<T, U>::setPlayer2Pos(Vec2D const pos)
{
    pl2.startPos = pos;
}

template<class T, class U>
PlayerState<T>const& MicroSimulation<T,U>::getPlayer1() const
{
    return pl1;
}

template<class T, class U>
PlayerState<U>const& MicroSimulation<T,U>::getPlayer2() const
{
    return pl2;
}


template<class T, class U>
bool MicroSimulation<T, U>::run(int const steps)
{
    for(int i = 0; i < steps; ++i)
    {
        timestep();
    }
    if(pl1.unitCount == 0 || pl2.unitCount == 0)
    {
        return true;
    }
    return false;
}

template <class T, class U>
void MicroSimulation<T, U>::setTimeSteps(size_t timeSteps)
{
    mTimeSteps = timeSteps;
}

template <class T, class U>
void MicroSimulation<T, U>::setTimeSlice(int timeSlice)
{
    mTimeSlice = timeSlice;
    for(auto unit : pl1.unitList)
    {
        unit->setTimeSlice(timeSlice);
    }
    for(auto unit : pl2.unitList)
    {
        unit->setTimeSlice(timeSlice);
    }
}


template <class T, class U>
void MicroSimulation<T, U>::timestep()
{
    pl1.timestep(pl2);
    pl2.timestep(pl1);
}

template <class T, class U>
Fitness MicroSimulation<T, U>::run(bool const reset)
{

    for(int i = 0; i < mTimeSteps; ++i)
    {
        if(pl1.unitCount == 0 || pl2.unitCount == 0)
        {
            break;
        }
        timestep();

    }
    Fitness res;
    double maxHealth = 0;;
    double maxMinerals_alive = 0;
    double maxGas_alive = 0;
    for(auto unit : pl1.unitList)
    {
        maxHealth += unit->getMaxHealth() + unit->getMaxShield();
        maxMinerals_alive += unit->getMinerals();
        maxGas_alive += unit->getGas();
        res.health += unit->getHealth() + unit->getShield();
        if(!unit->isDead())
        {
            res.health_alive += unit->getMaxHealth() + unit->getShield();
            res.minerals_alive += unit->getMinerals();
            res.gas_alive += unit->getGas();
        }
    }
    double maxDamage = 0;
    double maxMinerals_killed = 0;
    double maxGas_killed = 0;
    for(auto unit : pl2.unitList)
    {
        maxDamage += unit->getMaxHealth() + unit->getMaxShield();
        maxMinerals_killed += unit->getMinerals();
        maxGas_killed += unit->getGas();
        res.damage += unit->getHealth() + unit->getShield() - unit->getHealth() - unit->getShield();
        if(unit->isDead())
        {
            res.damage_killed += unit->getMaxHealth() + unit->getShield();
            res.minerals_killed += unit->getMinerals();
            res.gas_killed += unit->getGas();
        }
    }

    double maxDamage_inv = 1.0/maxDamage;
    res.damage *= maxDamage_inv;
    res.damage_killed *= maxDamage_inv;
    res.minerals_killed /= maxMinerals_killed;
    res.gas_killed /= maxGas_killed;

    double maxHealth_inv = 1.0/maxHealth;
    res.health *= maxHealth_inv;
    res.health_alive *= maxHealth_inv;
    res.minerals_alive /= maxMinerals_alive;
    res.gas_alive /= maxGas_alive;

    res.score = res.damage * res.damage_killed + res.minerals_killed + res.gas_killed;
    res.score += res.health * res.health_alive + res.minerals_alive + res.gas_alive;
    res.score /= 6.0;
    res.score *= 100.0;

    if(reset)
    {
        resetBothPlayers();
    }

    return res;
}





template<class T, class U>
void MicroSimulation<T,U>::setTracking(bool const tracking)
{
    if(tracking != mTracking)
    {
        for(auto unit : pl1.unitList)
        {
            unit->setTracking(tracking);
        }
        for(auto unit : pl2.unitList)
        {
            unit->setTracking(tracking);
        }
    }
}

template<class T, class U>
void MicroSimulation<T,U>::setTracking(bool const tracking, size_t const steps)
{
    setTracking(tracking);
    if(tracking)
    {
        for(auto unit : pl1.unitList)
        {
            unit->reservePathStorage(steps);
        }
        for(auto unit : pl2.unitList)
        {
            unit->reservePathStorage(steps);
        }
    }

}

template<class T, class U>
void MicroSimulation<T,U>::clearUnitPaths()
{
    for(auto unit : pl1.unitList)
    {
        unit->clearPath();
    }
    for(auto unit : pl2.unitList)
    {
        unit->clearPath();
    }
}

template<class T, class U>
size_t MicroSimulation<T,U>::getPlayer1ChromosomeLength() const
{
    return pl1.NGenes;
}

template<class T, class U>
size_t MicroSimulation<T,U>::getPlayer2ChromosomeLength() const
{
    return pl2.NGenes;
}

