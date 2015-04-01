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
MicroSimulation<T, U>::MicroSimulation(const pair<double,double> minPos, const pair<double,double> maxPos, const string& filePath1, const string& filePath2)
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
pair<double, double> MicroSimulation<T, U>::getMinPos() const
{
    return mMinPos;
}
template <class T, class U>
pair<double, double> MicroSimulation<T, U>::getMaxPos() const
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
void MicroSimulation<T,U>::initPlayer1(vector<string>const& unitList, UnitGenes const& genes, std::function<pair<double,double>(BaseUnit const&, BaseUnit const&)> friendForce, std::function<pair<double,double>(BaseUnit const&,BaseUnit const&)> enemyForce)
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

    auto pl1FuncMinX = [&](pair<double,double> const& pos,typename T::BUT const& unit)
        {
            double const dist = fabs(pos.first-unit.getX());
            if(dist < unit.getSpeed())
            {
                return pair<double,double>(LIMIT,0);
            }
            double const z = 1/pow(dist,3)*LIMIT;//+limit;

            return pair<double,double>(z,0);
        };
    auto pl1FuncMinY = [&](pair<double,double> const& pos,typename T::BUT const& unit)
        {
            double const dist = fabs(pos.second-unit.getY());
            if(dist < unit.getSpeed())
            {
                return pair<double,double>(0,LIMIT);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return pair<double,double>(0,z);
        };
    auto pl1FuncMaxX = [&](pair<double,double> const& pos,typename T::BUT const& unit)
        {
            double const dist = fabs(pos.first-unit.getX());
            if(dist < unit.getSpeed())
            {
                return pair<double,double>(-LIMIT,0);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return pair<double,double>(-z,0);
        };
    auto pl1FuncMaxY = [&](pair<double,double> const& pos,typename T::BUT const& unit)
        {
            double const dist = fabs(pos.second-unit.getY());
            if(dist < unit.getSpeed())
            {
                return pair<double,double>(0,-LIMIT);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return pair<double,double>(0,-z);
        };

    auto pl2FuncMinX = [&](pair<double,double> const& pos,typename U::BUT const& unit)
        {
            double dist = fabs(pos.first-unit.getX());
            if(dist < unit.getSpeed())
            {
                return pair<double,double>(LIMIT,0);
            }
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            //double const z = -pow(dist,4)+limit;
            return pair<double,double>(z,0);
        };
    auto pl2FuncMinY = [&](pair<double,double> const& pos,typename U::BUT const& unit)
        {
            double const dist = fabs(pos.second-unit.getY());
            if(dist < unit.getSpeed())
            {
                return pair<double,double>(0,LIMIT);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return pair<double,double>(0,z);
        };
    auto pl2FuncMaxX = [&](pair<double,double> const& pos,typename U::BUT const& unit)
        {
            double const dist = fabs(pos.first-unit.getX());
            if(dist < unit.getSpeed())
            {
                return pair<double,double>(-LIMIT,0);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return pair<double,double>(-z,0);
        };
    auto pl2FuncMaxY = [&](pair<double,double> const& pos,typename U::BUT const& unit)
        {
            double const dist = fabs(pos.second-unit.getY());
            if(dist < unit.getSpeed())
            {
                return pair<double,double>(0,-LIMIT);
            }
            //double const z = -pow(dist,4)+limit;
            double const z = 1/pow(dist,3)*LIMIT;//+limit;
            return pair<double,double>(0,-z);
        };

    pl1.potentialList.emplace_back(pl1.minPos,pl1FuncMinX);
    pl1.potentialList.emplace_back(pl1.maxPos,pl1FuncMaxX);
    pl1.potentialList.emplace_back(pl1.minPos,pl1FuncMinY);
    pl1.potentialList.emplace_back(pl1.maxPos,pl1FuncMaxY);
    pl2.potentialList.emplace_back(pl2.minPos,pl2FuncMinX);
    pl2.potentialList.emplace_back(pl2.maxPos,pl2FuncMaxX);
    pl2.potentialList.emplace_back(pl2.minPos,pl2FuncMinY);
    pl2.potentialList.emplace_back(pl2.maxPos,pl2FuncMaxY);
    double const pl1FieldSizeX = pl1.maxPos.first-pl1.minPos.first;
    double const pl1FieldSizeY = pl1.maxPos.second-pl1.minPos.second;
    double const pl2FieldSizeX = pl2.maxPos.first-pl2.minPos.first;
    double const pl2FieldSizeY = pl2.maxPos.second-pl2.minPos.second;
    pl1.startPos = pair<double,double>(pl1.minPos.first+pl1FieldSizeX/20.,pl1.minPos.second);
    pl2.startPos = pair<double,double>(pl2.maxPos.first-pl2FieldSizeX/20.,pl2.maxPos.second);

    for(auto unit : pl1.unitList)
    {
        pl1.startPos.second += pl1FieldSizeY/pl1.unitList.size();
        unit->setPos(pl1.startPos);
    }
    for(auto unit : pl2.unitList)
    {
        pl2.startPos.second -= pl2FieldSizeY/pl2.unitList.size();
        unit->setPos(pl2.startPos);
    }
}

template <class T, class U>
void MicroSimulation<T, U>::initBothPlayers(const vector<string>& unitList1, const vector<string>& unitList2)
{
    init1.init(unitList1, pl1);
    init2.init(unitList2, pl2);
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

template<class T, class U>
template<class UnitType>
double MicroSimulation<T, U>::sumUnitListCosts(const list<UnitType>& units)
{

    double sum = 0.;
    for (const UnitType& unit : units)
    {
        sum += unit.getMinerals();
        sum += GASTOMINERALS*unit.getGas();
    }
    return sum;
}

template<class T, class U>
double MicroSimulation<T, U>::sumPlayer1UnitCosts()
{
    double sum = 0;
    for(auto unitPtr : pl1.unitList)
    {
        sum += unitPtr->getMinerals();
        sum += GASTOMINERALS*unitPtr->getGas();
    }
    return sum;
}

template<class T, class U>
double MicroSimulation<T, U>::sumPlayer2UnitCosts()
{
    double sum = 0;
    for(auto unitPtr : pl2.unitList)
    {
        sum += unitPtr->getMinerals();
        sum += GASTOMINERALS*unitPtr->getGas();
    }
    return sum;
}

template <class T, class U>
double MicroSimulation<T, U>::compareCosts()
{
    double sum1 = sumPlayer1UnitCosts();
    double sum2 = sumPlayer2UnitCosts();
    return sum1 - sum2;
}


template<class T, class U>
PlayerStats MicroSimulation<T, U>::sumPlayer1UnitStats()
{
    PlayerStats stats;
    for(auto unitPtr : pl1.unitList)
    {
        if(unitPtr->isAirUnit())
        {
            stats.airArmor += unitPtr->getArmor();
            stats.airHealth += unitPtr->getHealth();
        }
        else
        {
            stats.groundArmor += unitPtr->getArmor();
            stats.groundHealth += unitPtr->getHealth();
        }
        stats.airDps += unitPtr->getAdps();
        stats.groundDps += unitPtr->getGdps();
    }
    return stats;
}

template<class T, class U>
PlayerStats MicroSimulation<T, U>::sumPlayer2UnitStats()
{
    PlayerStats stats;
    for(auto unitPtr : pl2.unitList)
    {
        if(unitPtr->isAirUnit())
        {
            stats.airArmor += unitPtr->getArmor();
            stats.airHealth += unitPtr->getHealth();
        }
        else
        {
            stats.groundArmor += unitPtr->getArmor();
            stats.groundHealth += unitPtr->getHealth();
        }
        stats.airDps += unitPtr->getAdps();
        stats.groundDps += unitPtr->getGdps();
    }
    return stats;
}


template <class T, class U>
int MicroSimulation<T, U>::compareStats()
{
    PlayerStats stats1 = sumPlayer1UnitStats();
    PlayerStats stats2 = sumPlayer2UnitStats();
    double oldGround1, oldGround2, oldAir1, oldAir2, subGround1, subGround2, subAir1, subAir2;
    for(;;)
    {
        if(stats1.groundHealth <= 0 && stats2.groundHealth <= 0 && stats1.groundHealth <= 0 && stats2.airHealth <= 0)
        {
            return 0;
        }
        else if(stats1.groundHealth <= 0 && stats1.airHealth <= 0)
        {
            return -1;
        }
        else if(stats2.groundHealth <= 0 && stats2.airHealth <= 0)
        {
            return 1;
        }
        oldGround1 = stats1.groundHealth;
        oldAir1 = stats1.airHealth;
        oldGround2 = stats2.groundHealth;
        oldAir2 = stats2.airHealth;
        subGround1 = stats2.groundDps-stats1.groundArmor;
        subAir1 = stats2.airDps-stats2.airArmor;
        subGround2 = stats1.groundDps-stats2.groundArmor;
        subAir2 = stats1.airDps-stats1.airArmor;
        stats1.groundHealth = stats1.groundHealth-subGround1 > 0 ? stats1.groundHealth-subGround1 : 0;
        stats1.airHealth = stats1.airHealth-subAir1 > 0 ? stats1.airHealth-subAir1 : 0;
        stats2.groundHealth = stats2.groundHealth-subGround2 > 0 ? stats2.groundHealth-subGround2 : 0;
        stats2.airHealth = stats2.airHealth-subAir2 > 0 ? stats1.airHealth-subAir2 : 0;
        if(oldGround1 == stats1.groundHealth && oldGround2 == stats2.groundHealth && oldAir1 == stats1.airHealth && oldAir2 == stats2.airHealth)
        {
            return 0;
        }
    }
}

template <class T, class U>
void MicroSimulation<T, U>::setPlayer1Pos(std::pair<double,double> const pos)
{
    pl1.startPos = pos;
}

template <class T, class U>
void MicroSimulation<T, U>::setPlayer2Pos(std::pair<double,double> const pos)
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
bool MicroSimulation<T, U>::run(size_t const steps)
{
    bool pl1Finished = true, pl2Finished = true;
    for(size_t i = 0; i < steps; i += mTimeFrame)
    {
        pl1Finished = true;
        pl2Finished = true;
        for(auto unit : pl1.unitList)
        {
            unit->attack(pl2);
        }
        for(auto unit : pl2.unitList)
        {
            unit->attack(pl1);
        }
        for(auto unit : pl1.unitList)
        {
            if(unit->getHealth() > EPS)
            {
                pl1Finished = false;
            }
            else
            {
                continue;
            }
            unit->move(pl1,pl2);
            unit->regenerate();
        }
        for(auto& unit : pl2.unitList)
        {
            if(unit->getHealth() > EPS)
            {
                pl2Finished = false;
            }
            else
            {
                continue;
            }
            unit->move(pl2,pl1);
            unit->regenerate();
        }
    }
    collectGarbage();
    if(pl1.unitList.empty() || pl2.unitList.empty())
    {
        return true;
    }

    if(pl1Finished || pl2Finished)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <class T, class U>
void MicroSimulation<T, U>::setTimeSteps(size_t timeSteps)
{
    mTimeSteps = timeSteps;
}

template <class T, class U>
void MicroSimulation<T, U>::setTimeFrame(size_t timeFrame)
{
    mTimeFrame = timeFrame;
}

template <class T, class U>
Fitness MicroSimulation<T, U>::run(bool const reset)
{
    initPotentialFields();
    pl1.unitCount = pl1.unitList.size();
    pl2.unitCount = pl2.unitList.size();
    for(size_t i = 0; i < mTimeSteps; i += mTimeFrame)
    {
        if(pl1.unitCount == 0 || pl2.unitCount == 0)
        {
            break;
        }
        for(auto unit : pl1.unitList)
        {
            unit->attack(pl2);
        }
        for(auto unit : pl2.unitList)
        {
            unit->attack(pl1);
        }
        for(auto unit : pl1.unitList)
        {
            if(unit->getHealth() < EPS)
            {
                continue;
            }
            unit->move(pl1,pl2);
        }
        for(auto& unit : pl2.unitList)
        {
            if(unit->getHealth() < EPS)
            {
                continue;
            }
            unit->move(pl2,pl1);
        }
        pl1.regenerate();
        pl2.regenerate();
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
