#ifndef _MICROSIMULATION_H_
#define _MICROSIMULATION_H_

#include<memory>
#include<string>
#include<vector>
#include<list>
#include<utility>
#include<functional>

#include "PlayerState.h"
#include "Unit.h"
#include "InitPlayerUnits.h"
#include "Utilities.h"

using std::shared_ptr;
using std::string;
using std::vector;
using std::list;
using std::pair;

struct PlayerStats
{
PlayerStats()
: groundArmor(0), groundHealth(0), groundDps(0), airArmor(0), airHealth(0), airDps(0)
        {}

    double groundArmor;
    double groundHealth;
    double groundDps;
    double airArmor;
    double airHealth;
    double airDps;

};

template<class T, class U>
    class MicroSimulation
{
private:

    Vec2D mMinPos;
    Vec2D mMaxPos;
    string mFilePath1;
    string mFilePath2;
    InitPlayerUnits<T> init1;
    InitPlayerUnits<U> init2;
    PlayerState<T> pl1;
    PlayerState<U> pl2;
    bool mTracking = false;
    int mTimeSteps = 30000;
    int mTimeSlice = 10;
    template <class V> void collectPlayerGarbage(PlayerState<V>&);

public:
    MicroSimulation(MicroSimulation const& microSim);
    MicroSimulation(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2);
    void initPlayer1(vector<string> const&);
    void initPlayer1(vector<string>const&, UnitGenes const& genes, std::function<Vec2D(BaseUnit &, BaseUnit &)> friendFunc, std::function<Vec2D(BaseUnit &,BaseUnit &)> enemyFunc);
    void initPlayer2(vector<string> const&);
    void setPlayer1Genes(UnitGenes const& genes);
    void setPlayer2Genes(UnitGenes const& genes);
    void setPlayer1Pos(Vec2D const pos);
    void setPlayer2Pos(Vec2D const pos);
    void initBothPlayers(vector<string> const&, vector<string> const&);
    void clearPlayer1();
    void clearPlayer2();
    void clearBothPlayers();
    void resetPlayer1();
    void resetPlayer2();
    void resetBothPlayers();
    PlayerState<T>const& getPlayer1() const;
    PlayerState<U>const& getPlayer2() const;
    Vec2D getMinPos() const;
    Vec2D getMaxPos() const;
    string getFilePath1() const;
    string getFilePath2() const;
    bool run(int const steps);
    Fitness run(bool const reset);
    void setTimeSteps(size_t timeSteps);
    void setTimeSlice(int timeSlice);

    void timestep();

    string determineWinner();


    void collectGarbage();


    void setTracking(bool const tracking);
    void setTracking(bool const tracking, size_t const steps);
    void clearUnitPaths();

};

#endif
