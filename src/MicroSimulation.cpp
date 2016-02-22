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
#include<chrono>


template <class T, class U>
MicroSimulation<T, U>::MicroSimulation(const Vec2D minPos, const Vec2D maxPos, const std::string& filePath1, const std::string& filePath2)
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
std::string MicroSimulation<T, U>::getFilePath1() const
{
    return mInfoDirName1;
}

template <class T, class U>
std::string MicroSimulation<T, U>::getFilePath2() const
{
    return mInfoDirName2;
}



template <class T, class U>
void MicroSimulation<T, U>::initPlayer1(const std::vector<std::string>& unitList)
{
    init1.init(unitList, mInfoDirName1, pl1);
    setUnitStartPositions(pl1, pl1.minPos.x + 15);
}

template <class T, class U>
void MicroSimulation<T, U>::initPlayer2(const std::vector<std::string>& unitList)
{
    init2.init(unitList, mInfoDirName2, pl2);
    setUnitStartPositions(pl2, pl2.maxPos.x - 15);

}

template<class T, class U>
template<class W>
void MicroSimulation<T, U>::setUnitStartPositions(PlayerState<W>& pl, double const x_start)
{

    double const fieldSizeY = pl.maxPos.y-pl.minPos.y;
    double const x = x_start;
    double y1 = pl.minPos.y + 0.5*fieldSizeY;
    double y2 = y1;

    int i = 0;

    for(auto unit : pl.unitList)
    {
        if(i % 2 == 0)
        {
            y1 += 3.0*unit->getSize();
            unit->setStartPos(Vec2D(x, y1));
            y1 += 3.0*unit->getSize();
        }
        else
        {
            y2 -= 3.0*unit->getSize();
            unit->setStartPos(Vec2D(x, y2));
            y2 -= 3.0*unit->getSize();
        }
        unit->resetPos();
        ++i;
    }
}



template <class T, class U>
void MicroSimulation<T, U>::initBothPlayers(const std::vector<std::string>& unitList1, const std::vector<std::string>& unitList2)
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
unsigned long MicroSimulation<T, U>::getNumberOfRuns() const
{
    return mRuns;
}

template <class T, class U>
Fitness MicroSimulation<T, U>::run(bool const reset, Player const player)
{

    ++mRuns;
/*
    if(mTracking)
    {
        // Save the unit names
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

        mFile1.setf(std::ios::fixed,std::ios::floatfield);
        mFile1.precision(2);
        mFile2.setf(std::ios::fixed,std::ios::floatfield);
        mFile2.precision(2);
    }
*/
    for(int i = 0; i < mTimeSteps; ++i)
    {
        // If all units of one player are dead, stop the simulation
        if(pl1.unitCount == 0 || pl2.unitCount == 0)
        {
            break;
        }
        timestep();
        /*
        if(mTracking && mTimeSteps % 2 == 0)
        {
            // Save paths and sum of health and shield
            for(auto const unitPtr : pl1.unitList)
            {
                if(unitPtr->isDead())
                {
                    mFile1 << "-" << "\t\t\t";
                }
                else
                {
                    mFile1 << unitPtr->getX() << "," << unitPtr->getY() << "," << (unitPtr->getHealth()+unitPtr->getShield())/(unitPtr->getMaxHealth() + unitPtr->getMaxShield()) << "\t\t\t";
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
                    mFile2 << unitPtr->getX() << "," << unitPtr->getY() << "," << (unitPtr->getHealth()+unitPtr->getShield())/(unitPtr->getMaxHealth() + unitPtr->getMaxShield()) << "\t\t\t";
                }
            }
            mFile2 << std::endl;
        }*/
    }
    /*
    if(mTracking)
    {
        mFile1.close();
        mFile2.close();
    }*/

    // Evaluate the outcome
    Fitness res;
    double maxDamage = 0;
    double maxHealth = 0;

    if(player == Player::first)
    {
        for(auto unit : pl1.unitList)
        {
            maxHealth += (unit->getMaxHealth() + unit->getMaxShield())*(unit->getMinerals()+unit->getGas());
            if(unit->hasAttacked()) res.health += (unit->getHealth() + unit->getShield())*(unit->getMinerals()+unit->getGas());
        }

        for(auto unit : pl2.unitList)
        {
            maxDamage += (unit->getMaxHealth() + unit->getMaxShield())*(unit->getMinerals()+unit->getGas());
            res.damage += (unit->getMaxHealth() + unit->getMaxShield() - std::max(0.0, unit->getHealth()) - std::max(0.0, unit->getShield()))*(unit->getMinerals()+unit->getGas());;
        }

    }
    else
    {
        for(auto unit : pl2.unitList)
        {
            maxHealth += (unit->getMaxHealth() + unit->getMaxShield())*(unit->getMinerals()+unit->getGas());
            if(unit->hasAttacked()) res.health += (unit->getHealth() + unit->getShield())*(unit->getMinerals()+unit->getGas());
        }

        for(auto unit : pl1.unitList)
        {
            maxDamage += (unit->getMaxHealth() + unit->getMaxShield())*(unit->getMinerals()+unit->getGas());
            res.damage += (unit->getMaxHealth() + unit->getMaxShield() - std::max(0.0, unit->getHealth()) - std::max(0.0, unit->getShield()))*(unit->getMinerals()+unit->getGas());;
        }
    }

    res.damage /= maxDamage;
    res.health /= maxHealth;

    // Can be used as approximation for multiple objectives as single objective
    res.score = (res.damage + res.health) * 50.0;

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
void MicroSimulation<T,U>::enableTracking(std::string const& fileName1, std::string const& fileName2)
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



