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
using std::ios;

template <class T, class U>
MicroSimulation<T, U>::MicroSimulation(const Vec2D minPos, const Vec2D maxPos, const string& filePath1, const string& filePath2)
    : mMinPos(minPos), mMaxPos(maxPos), mInfoDirName1(filePath1), mInfoDirName2(filePath2), init1(filePath1), init2(filePath2)
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
    return mInfoDirName1;
}

template <class T, class U>
string MicroSimulation<T, U>::getFilePath2() const
{
    return mInfoDirName2;
}



template <class T, class U>
void MicroSimulation<T, U>::initPlayer1(const vector<string>& unitList)
{
    init1.init(unitList, mInfoDirName1, pl1);
    double const fieldSizeX = pl1.maxPos.x-pl1.minPos.x;
    double const fieldSizeY = pl1.maxPos.y-pl1.minPos.y;
    Vec2D startPos = Vec2D(pl1.minPos.x + 10, pl1.minPos.y + 0.5*fieldSizeY);
    int i = 0;
    for(auto unit : pl1.unitList)
    {
        if(i % 2 == 1)
        {
            unit->setStartPos(Vec2D(startPos.x, startPos.y + i*5));
        }
        else
        {
            unit->setStartPos(Vec2D(startPos.x, startPos.y - i*5));
        }
        unit->resetPos();
        ++i;
    }
}

template <class T, class U>
void MicroSimulation<T, U>::initPlayer2(const vector<string>& unitList)
{
    init2.init(unitList, mInfoDirName2, pl2);
    double const fieldSizeX = pl2.maxPos.x-pl2.minPos.x;
    double const fieldSizeY = pl2.maxPos.y-pl2.minPos.y;
    Vec2D startPos = Vec2D(pl2.maxPos.x - 10, pl2.minPos.y + 0.5*fieldSizeY);
    int i = 0;
    for(auto unit : pl2.unitList)
    {
        if(i % 2 == 1)
        {
            unit->setStartPos(Vec2D(startPos.x, startPos.y + i*5));
        }
        else
        {
            unit->setStartPos(Vec2D(startPos.x, startPos.y - i*5));
        }
        unit->resetPos();
        ++i;
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
Fitness MicroSimulation<T, U>::run(bool const reset, Player const player)
{

    if(mTracking)
    {
        mFile1.open(mTrackingFileName1);
        for(auto const unitPtr : pl1.unitList)
        {
            mFile1 << unitPtr->getName() << "\t\t\t";
        }
        mFile1 << std::endl;
        mFile2.open(mTrackingFileName2);
        for(auto const unitPtr : pl2.unitList)
        {
            mFile2 << unitPtr->getName() << "\t\t\t";
        }
        mFile2 << std::endl;

        for(auto const unitPtr : pl1.unitList)
        {
            mFile1 << unitPtr->getSize() << "\t\t\t";
        }
        mFile1 << std::endl;
        for(auto const unitPtr : pl2.unitList)
        {
            mFile2 << unitPtr->getSize() << "\t\t\t";
        }
        mFile2 << std::endl;

        mFile1.setf(ios::fixed,ios::floatfield);
        mFile1.precision(2);
        mFile2.setf(ios::fixed,ios::floatfield);
        mFile2.precision(2);
    }

    for(int i = 0; i < mTimeSteps; ++i)
    {
        if(pl1.unitCount == 0 || pl2.unitCount == 0)
        {
            break;
        }
        timestep();
        if(mTracking)
        {
            for(auto const unitPtr : pl1.unitList)
            {
                if(unitPtr->isDead())
                {
                    mFile1 << "-" << "\t\t\t";
                }
                else
                {
                    mFile1 << unitPtr->getX() << "," << unitPtr->getY() << "\t\t\t";
                }
            }
            mFile1 << std::endl;

            for(auto const unitPtr : pl2.unitList)
            {
                if(unitPtr->isDead())
                {
                    mFile2 << "-" << "\t\t\t";
                }
                else
                {
                    mFile2 << unitPtr->getX() << "," << unitPtr->getY() << "\t\t\t";
                }
            }
            mFile2 << std::endl;
        }
    }
    if(mTracking)
    {
        mFile1.close();
        mFile2.close();
    }

    Fitness res;
    double maxDamage = 0;
    double maxHealth = 0;

    if(player == Player::first)
    {
        for(auto unit : pl1.unitList)
        {
            maxHealth += (unit->getMaxHealth() + unit->getMaxShield())*(unit->getMinerals()+unit->getGas());
            res.health += (unit->getHealth() + unit->getShield())*(unit->getMinerals()+unit->getGas());
        }

        for(auto unit : pl2.unitList)
        {
            maxDamage += (unit->getMaxHealth() + unit->getMaxShield())*(unit->getMinerals()+unit->getGas());
            res.damage += (unit->getMaxHealth() + unit->getMaxShield() - unit->getHealth() - unit->getShield())*(unit->getMinerals()+unit->getGas());;
        }


    }
    else
    {
        for(auto unit : pl2.unitList)
        {
            maxHealth += (unit->getMaxHealth() + unit->getMaxShield())*(unit->getMinerals()+unit->getGas());
            res.health += (unit->getHealth() + unit->getShield())*(unit->getMinerals()+unit->getGas());
        }

        for(auto unit : pl1.unitList)
        {
            maxDamage += (unit->getMaxHealth() + unit->getMaxShield())*(unit->getMinerals()+unit->getGas());
            res.damage += (unit->getMaxHealth() + unit->getMaxShield() - unit->getHealth() - unit->getShield())*(unit->getMinerals()+unit->getGas());;
        }
    }

    res.damage /= maxDamage;

    res.health /= maxHealth;

    res.score = (res.damage + res.health) * 100.0 / 2.0;

    res.timeSteps = std::max(pl1.time / pl1.timeSlice, pl2.time / pl2.timeSlice);

    if(reset)
    {
        resetBothPlayers();
    }

    return res;
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

template<class T, class U>
void MicroSimulation<T,U>::enableTracking(string const& fileName1, string const& fileName2)
{
    mTracking = true;
    mTrackingFileName1 = fileName1;
    mTrackingFileName2 = fileName2;
}

template<class T, class U>
void MicroSimulation<T,U>::disableTracking()
{
    mTracking = false;
}



