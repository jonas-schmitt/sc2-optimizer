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
    init1.init(unitList, pl1);
    pl1.unitCount = pl1.unitList.size();
}
template <class T, class U>
void MicroSimulation<T,U>::initPlayer1(vector<string>const& unitList, UnitGenes const& genes, std::function<Vec2D(BaseUnit &, BaseUnit &)> friendForce, std::function<Vec2D(BaseUnit &,BaseUnit &)> enemyForce)
{
    init1.init(unitList,pl1);
    for(auto& unit : pl1.unitList)
    {
        unit->setGenes(genes);
        unit->setFriendForce(friendForce);
        unit->setEnemyForce(enemyForce);
    }
    pl1.unitCount = pl1.unitList.size();
}

template <class T, class U>
void MicroSimulation<T, U>::initPlayer2(const vector<string>& unitList)
{
    init2.init(unitList, pl2);
    pl2.unitCount = pl2.unitList.size();
}

template <class T, class U>
void MicroSimulation<T, U>::initPotentialFields()
{

    auto pl1FuncMinX = [&](Vec2D const& pos,typename T::BUT const& unit)
        {
            double const dist = fabs(pos.x-unit.getX());
            if(dist < unit.getSpeed())
            {
                return Vec2D(LIMIT,0);
            }
            double const z = 1/pow(dist,3)*LIMIT;//+limit;

            return Vec2D(z,0);
        };
    auto pl1FuncMinY = [&](Vec2D const& pos,typename T::BUT const& unit)
        {
            double const dist = fabs(pos.y-unit.getY());
            if(dist < unit.getSpeed())
            {
                return Vec2D(0,LIMIT);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return Vec2D(0,z);
        };
    auto pl1FuncMaxX = [&](Vec2D const& pos,typename T::BUT const& unit)
        {
            double const dist = fabs(pos.x-unit.getX());
            if(dist < unit.getSpeed())
            {
                return Vec2D(-LIMIT,0);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return Vec2D(-z,0);
        };
    auto pl1FuncMaxY = [&](Vec2D const& pos,typename T::BUT const& unit)
        {
            double const dist = fabs(pos.y-unit.getY());
            if(dist < unit.getSpeed())
            {
                return Vec2D(0,-LIMIT);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return Vec2D(0,-z);
        };

    auto pl2FuncMinX = [&](Vec2D const& pos,typename U::BUT const& unit)
        {
            double dist = fabs(pos.x-unit.getX());
            if(dist < unit.getSpeed())
            {
                return Vec2D(LIMIT,0);
            }
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            //double const z = -pow(dist,4)+limit;
            return Vec2D(z,0);
        };
    auto pl2FuncMinY = [&](Vec2D const& pos,typename U::BUT const& unit)
        {
            double const dist = fabs(pos.y-unit.getY());
            if(dist < unit.getSpeed())
            {
                return Vec2D(0,LIMIT);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return Vec2D(0,z);
        };
    auto pl2FuncMaxX = [&](Vec2D const& pos,typename U::BUT const& unit)
        {
            double const dist = fabs(pos.x-unit.getX());
            if(dist < unit.getSpeed())
            {
                return Vec2D(-LIMIT,0);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return Vec2D(-z,0);
        };
    auto pl2FuncMaxY = [&](Vec2D const& pos,typename U::BUT const& unit)
        {
            double const dist = fabs(pos.y-unit.getY());
            if(dist < unit.getSpeed())
            {
                return Vec2D(0,-LIMIT);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return Vec2D(0,-z);
        };

    pl1.potentialList.emplace_back(pl1.minPos,pl1FuncMinX);
    pl1.potentialList.emplace_back(pl1.maxPos,pl1FuncMaxX);
    pl1.potentialList.emplace_back(pl1.minPos,pl1FuncMinY);
    pl1.potentialList.emplace_back(pl1.maxPos,pl1FuncMaxY);
    pl2.potentialList.emplace_back(pl2.minPos,pl2FuncMinX);
    pl2.potentialList.emplace_back(pl2.maxPos,pl2FuncMaxX);
    pl2.potentialList.emplace_back(pl2.minPos,pl2FuncMinY);
    pl2.potentialList.emplace_back(pl2.maxPos,pl2FuncMaxY);
    double const pl1FieldSizeX = pl1.maxPos.x-pl1.minPos.x;
    double const pl1FieldSizeY = pl1.maxPos.y-pl1.minPos.y;
    double const pl2FieldSizeX = pl2.maxPos.x-pl2.minPos.x;
    double const pl2FieldSizeY = pl2.maxPos.y-pl2.minPos.y;
    pl1.startPos = Vec2D(pl1.minPos.x+pl1FieldSizeX/20.,pl1.minPos.y);
    pl2.startPos = Vec2D(pl2.maxPos.x-pl2FieldSizeX/20.,pl2.maxPos.y);

    for(auto unit : pl1.unitList)
    {
        pl1.startPos.y += pl1FieldSizeY/pl1.unitList.size();
        unit->setPos(pl1.startPos);
    }
    for(auto unit : pl2.unitList)
    {
        pl2.startPos.y -= pl2FieldSizeY/pl2.unitList.size();
        unit->setPos(pl2.startPos);
    }
}

template <class T, class U>
void MicroSimulation<T, U>::initBothPlayers(const vector<string>& unitList1, const vector<string>& unitList2)
{
    init1.init(unitList1, pl1);
    init2.init(unitList2, pl2);
    resetTime();
    pl1.unitCount = pl1.unitList.size();
    pl2.unitCount = pl2.unitList.size();
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
void MicroSimulation<T, U>::resetTime()
{
    pl1.time = 0;
    pl1.movementTimer= pl1.movementUpdate;
    pl1.regenerationTimer = pl1.regenerationUpdate;
    pl2.time = 0;
    pl2.movementTimer= pl2.movementUpdate;
    pl2.regenerationTimer = pl2.regenerationUpdate;
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
    initPotentialFields();
    resetTime();
    pl1.unitCount = pl1.unitList.size();
    pl2.unitCount = pl2.unitList.size();

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
    if(res.damage < EPS)
    {
        res.score = 0;
    }
    else
    {
        res.score = res.health+res.damage;
    }
    res.health /= maxHealth;
    res.damage /= maxDamage;
    return res;
}

template <class T, class U>
void MicroSimulation<T,U>::collectGarbage()
{

//    auto tmp = std::async(std::launch::async, [&] ()
//    {
    collectPlayerGarbage(pl1);
//    });
    collectPlayerGarbage(pl2);
//    tmp.get();
}

template <class T, class U>
template <class V>
void MicroSimulation<T,U>::collectPlayerGarbage(PlayerState<V>& pl)
{
    //    auto tmp1 = std::async(std::launch::async,
    //        [&] ()
    //    {
    pl.removeZombies(pl.unitList0);
    pl.removeZombies(pl.unitList1);
    pl.removeZombies(pl.unitList2);
    pl.removeZombies(pl.unitList3);
    pl.removeZombies(pl.unitList4);
    pl.removeZombies(pl.unitList5);
    //    });

    //    auto tmp2 = std::async(std::launch::async,
    //        [&] ()
    //    {
    pl.removeZombies(pl.unitList6);
    pl.removeZombies(pl.unitList7);
    pl.removeZombies(pl.unitList8);
    pl.removeZombies(pl.unitList9);
    pl.removeZombies(pl.unitList10);
    //    });

    pl.removeZombies(pl.unitList11);
    pl.removeZombies(pl.unitList12);
    pl.removeZombies(pl.unitList13);
    pl.removeZombies(pl.unitList14);
    pl.removeZombies(pl.unitList15);
    pl.removeZombies(pl.unitList16);
    pl.removeZombies(pl.unitList17);
    //    tmp1.get();
    //    tmp2.get();
}

template<class T, class U>
string MicroSimulation<T,U>::determineWinner()
{
    collectGarbage();
    if(pl1.unitList.empty())
    {
        return "Player2";
    }
    else if(pl2.unitList.empty())
    {
        return "Player1";
    }
    else
    {
        return "Nobody";
    }
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
