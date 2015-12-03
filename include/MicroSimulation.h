#ifndef _MICROSIMULATION_H_
#define _MICROSIMULATION_H_

#include<memory>
#include<string>
#include<vector>

#include<utility>
#include<functional>
#include<fstream>

#include "PlayerState.h"
#include "Unit.h"
#include "InitPlayerUnits.h"
#include "Utilities.h"

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


// Class representing the simulation of an encounter
template<class T, class U>
    class MicroSimulation final
{
private:

    Vec2D mMinPos, mMaxPos;

    // directories containing the information of the units of both players
    std::string mInfoDirName1, mInfoDirName2;

    InitPlayerUnits<T> init1;
    InitPlayerUnits<U> init2;
    PlayerState<T> pl1;
    PlayerState<U> pl2;
    bool mTracking = false;
    // Maximum number of time steps
    long mTimeSteps = 12000;
    // Length of a time step
    int mTimeSlice = 10;

    // Number of executed runs
    unsigned long mRuns = 0;

    // Names of the files where the unit paths should be stored
    std::string mTrackingFileName1, mTrackingFileName2;
    // Corresponding filestreams
    std::ofstream mFile1, mFile2;



public:
    // Copy constructor
    MicroSimulation(MicroSimulation const& microSim);

    // Main constructor
    // filePath1, filePath2: Paths to the directories containing the unit stats
    MicroSimulation(Vec2D const minPos, Vec2D const maxPos, std::string const& filePath1, std::string const& filePath2);
    template<typename V> void setPlayerChromosome(PlayerState<V>& pl, Chromosome const&);
    void setPlayer1Chromosome(Chromosome const&);
    void setPlayer2Chromosome(Chromosome const&);
    void setPlayer1Pos(Vec2D const pos);
    void setPlayer2Pos(Vec2D const pos);

    void initPlayer1(std::vector<std::string> const&);
    void initPlayer2(std::vector<std::string> const&);
    void initBothPlayers(std::vector<std::string> const&, std::vector<std::string> const&);

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
    std::string getFilePath1() const;
    std::string getFilePath2() const;

    // Main function for running a simulation
    // Player denotes which player should be evaluated (either Player::Player1 or
    Fitness run(bool const reset, Player const player);

    void setTimeSteps(size_t timeSteps);
    void setTimeSlice(int timeSlice);

    // Time stepping method used internally within the run method
    void timestep();

    // Predict winner based on the unit stats
    std::string determineWinner();



    // Enable or disable the tracking of unit paths
    void enableTracking(std::string const& fileName1, std::string const& fileName2);
    void disableTracking();
    void clearUnitPaths();

    size_t getPlayer1ChromosomeLength() const;
    size_t getPlayer2ChromosomeLength() const;

    template<class W>
    void setUnitStartPositions(PlayerState<W>& pl, double const x_start);
    unsigned long getNumberOfRuns() const;
};

#endif
