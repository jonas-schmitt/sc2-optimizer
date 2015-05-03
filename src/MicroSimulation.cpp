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
void MicroSimulation<T,U>::initPlayer1(vector<string>const& unitList, UnitGenes const& genes, std::function<Vec2D(BaseUnit &, BaseUnit &)> friendForce, std::function<Vec2D(BaseUnit &,BaseUnit &)> enemyForce)
{
    initPlayer1(unitList);
    for(auto& unit : pl1.unitList)
    {
        unit->setGenes(genes);
        unit->setFriendForce(friendForce);
        unit->setEnemyForce(enemyForce);
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
void MicroSimulation<T,U>::setPlayer1Genes(UnitGenes const& genes)
{
    for(auto unit : pl1.unitList)
    {
        unit->setGenes(genes);
    }
}

template<class T, class U>
void MicroSimulation<T,U>::setPlayer2Genes(UnitGenes const& genes)
{
    for(auto unit : pl2.unitList)
    {
        unit->setGenes(genes);
    }
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
    double maxHealth = 0;
    double maxDamage = 0;
    for(auto unit : pl1.unitList)
    {
        maxHealth += unit->getMaxHealth();
        maxHealth += unit->getMaxShield();
        res.health += unit->getHealth();
        res.health += unit->getShield();

    }
    for(auto unit : pl2.unitList)
    {
        maxDamage += unit->getMaxHealth();
        maxDamage += unit->getMaxShield();

        res.damage += (unit->getMaxHealth() - unit->getHealth());
        res.damage += (unit->getMaxShield() - unit->getShield());
    }

    if(reset)
    {
        resetBothPlayers();
    }

    res.score = res.health+res.damage;
    res.health /= maxHealth;
    res.damage /= maxDamage;
    if(res.damage < res.health)
    {
        res.score = 0;
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
