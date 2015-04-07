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
    pair<double,double> mMinPos;
    pair<double,double> mMaxPos;
    string mFilePath1;
    string mFilePath2;
    InitPlayerUnits<T> init1;
    InitPlayerUnits<U> init2;
    PlayerState<T> pl1;
    PlayerState<U> pl2;
    bool mTracking = false;
    int mTimeSteps = 10000;
    int mTimeSlice = 10;
    template <class V> void collectPlayerGarbage(PlayerState<V>&);

public:
    MicroSimulation(MicroSimulation const& microSim);
    MicroSimulation(pair<double, double> const minPos, pair<double, double> const maxPos, string const& filePath1, string const& filePath2);
    void initPlayer1(vector<string> const&);
    void initPlayer1(vector<string>const&, UnitGenes const& genes, std::function<pair<double,double>(BaseUnit const&, BaseUnit const&)> friendFunc, std::function<pair<double,double>(BaseUnit const&,BaseUnit const&)> enemyFunc);
    void initPlayer2(vector<string> const&);
    void setPlayer1Genes(UnitGenes const& genes);
    void setPlayer2Genes(UnitGenes const& genes);
    void setPlayer1Pos(pair<double,double> const pos);
    void setPlayer2Pos(pair<double,double> const pos);
    void initBothPlayers(vector<string> const&, vector<string> const&);
    void clearPlayer1();
    void clearPlayer2();
    void clearBothPlayers();
    void resetPlayer1();
    void resetPlayer2();
    void resetBothPlayers();
    PlayerState<T>const& getPlayer1() const;
    PlayerState<U>const& getPlayer2() const;
    pair<double,double> getMinPos() const;
    pair<double,double> getMaxPos() const;
    string getFilePath1() const;
    string getFilePath2() const;
    bool run(int const steps);
    Fitness run(bool const reset);
    void setTimeSteps(size_t timeSteps);
    void setTimeSlice(int timeSlice);

    void timestep();

    string determineWinner();

    // returns a positive value if player1's units cost more resources than player2's
    double compareCosts();

    // returns 1 if player1 wins, returns -1 if player two wins, returns 0 in case of a draw
    int compareStats();

    void collectGarbage();

    void initPotentialFields();

    void setTracking(bool const tracking);
    void setTracking(bool const tracking, size_t const steps);
    void clearUnitPaths();

};

#endif
