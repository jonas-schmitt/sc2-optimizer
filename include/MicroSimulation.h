#ifndef _MICROSIMULATION_H_
#define _MICROSIMULATION_H_

#include<memory>
#include<string>
#include<vector>
#include<list>
#include<utility>
#include<functional>
#include<fstream>

#include "PlayerState.h"
#include "Unit.h"
#include "InitPlayerUnits.h"
#include "Utilities.h"

using std::shared_ptr;
using std::string;
using std::vector;
using std::list;
using std::pair;
using std::ofstream;

struct PlayerStats final
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
    class MicroSimulation final
{
private:


    Vec2D mMinPos, mMaxPos;
    string mInfoDirName1, mInfoDirName2;
    InitPlayerUnits<T> init1;
    InitPlayerUnits<U> init2;
    PlayerState<T> pl1;
    PlayerState<U> pl2;
    bool mTracking = false;
    long mTimeSteps = 30000;
    int mTimeSlice = 10;

    string mTrackingFileName1, mTrackingFileName2;

    ofstream mFile1, mFile2;

public:
    MicroSimulation(MicroSimulation const& microSim);
    MicroSimulation(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2);
    void initPlayer1(vector<string> const&);
    void initPlayer2(vector<string> const&);
    template<typename V> void setPlayerChromosome(PlayerState<V>& pl, Chromosome const&);
    void setPlayer1Chromosome(Chromosome const&);
    void setPlayer2Chromosome(Chromosome const&);
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
    Fitness run(bool const reset, Player const player);
    void setTimeSteps(size_t timeSteps);
    void setTimeSlice(int timeSlice);

    void timestep();

    string determineWinner();



    void enableTracking(string const& fileName1, string const& fileName2);
    void disableTracking();

    void clearUnitPaths();

    size_t getPlayer1ChromosomeLength() const;
    size_t getPlayer2ChromosomeLength() const;

};

#endif
