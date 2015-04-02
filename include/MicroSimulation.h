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
    size_t mTimeSteps = 450000;
    int mTimeSlice = 10;
    template<typename UnitType> double sumUnitListCosts(list<UnitType> const& units);
    double sumPlayer1UnitCosts();
    double sumPlayer2UnitCosts();
    PlayerStats sumPlayer1UnitStats();
    PlayerStats sumPlayer2UnitStats();
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
    bool run(size_t steps);
    Fitness run(bool const reset);
    void setTimeSteps(size_t timeSteps);
    void setTimeSlice(int timeSlice);

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

    void writePotentialToOut(std::string out)
    {
        std::ofstream outF(out.c_str(), std::ofstream::out);
        auto in = getPlayer1().unitList.front();
        std::vector<std::vector<double>> a;
        for (size_t i = 0; i < 180; ++i)
        {
            std::vector<double> b;
            for (size_t j = 0; j < 180; ++j)
            {
                b.push_back(0);
            }
            a.push_back(b);
        }
        //for (auto &in : getPlayer1().unitList)
        std::vector<std::vector<double>> bla;
        for (size_t i = 0; i < 150; i = i+1)
        {
            std::vector<double> t;
            for (size_t j = 0; j < 150; ++j)
            {
                in->setPos(i,j);
                t.push_back(in->getPot(pl1,pl2,a));
            }
            bla.push_back(t);
        }
        //outF << unit->getPos().first << "\t" << unit->getPos().second << "\t" << unit->getPot(pl1, pl2);
        //outF << std::endl;
        std::cout << "fin" <<std::endl;
        for (size_t i = 0; i < 150; i = i+2)
        {
            for (size_t j = 0; j < 150; ++j)
            {
                outF << i << "\t" << j << "\t" << bla.at(j).at(i);
                outF << std::endl;
            }
            for (size_t j = 150; j > 0; --j)
            {
                outF << i << "\t" << j-1 << "\t" << bla.at(j-1).at(i);
                outF << std::endl;
            }
        }
        outF.close();
        /*
          for (size_t i = 0; i < 150; i = i+2)
          {
          for (size_t j = 0; j < 150; ++j)
          {
          unit->setPos(i,j);
          outF << i << "\t" << j << "\t" << unit->getPot(pl1, pl2);
          outF << std::endl;
          }
          for (size_t j = 150; j > 0; --j)
          {
          unit->setPos(i,j-1);
          outF << i << "\t" << j-1 << "\t" << unit->getPot(pl1, pl2);
          outF << std::endl;
          }
          }
        */
        outF.close();
    }

    void getPotentialField()
    {
        size_t i = 0;
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
        run((size_t)5);
        writePotentialToOut("./potentials/potential"+std::to_string(i++)+".dat");
    }

};

#endif
