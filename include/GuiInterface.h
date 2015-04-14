#ifndef _GUI_INTERFACE_H
#define _GUI_INTERFACE_H

#include <iostream>
#include <unistd.h>

#include "MicroSimulation.h"
#include "Unit.h"
#include "BaseSelector.h"

//class GuiWindow; //forward declaration

template<class T, class U>
    class GuiInterface : public BaseSelector, MicroSimulation<T,U>
{
private:
    GuiWindow *mGui;

public:
    GuiInterface(Vec2D minPos, Vec2D maxPos, const string& filePath1, const string& filePath2);
    void setGuiWindow(GuiWindow *gui);
    bool run(int const steps);
    void run();
    void reset();

    void initBothPlayers(const vector<string>&, const vector<string>&);

    PlayerState<T>const& getPlayer1() const;
    PlayerState<U>const& getPlayer2() const;

    SimulationResult getPlayer1Result();
    SimulationResult getPlayer2Result();

    void clearBothPlayers();
    void setPlayer1Genes(UnitGenes const& genes) {MicroSimulation<T,U>::setPlayer1Genes(genes);}
    void setPlayer2Genes(UnitGenes const& genes) {MicroSimulation<T,U>::setPlayer2Genes(genes);}

    void setPositions(std::pair<int, int> &, std::pair<int, int> &);

    void collectGarbage();
    void setTracking(bool const tracking);
    void setTracking(bool const tracking, size_t const steps);
    void setTimeSteps(size_t steps);
    void clearUnitPaths();


};

#endif
