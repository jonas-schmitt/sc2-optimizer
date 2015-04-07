#include "GuiInterface.h"

template<class T, class U>
GuiInterface<T, U>::GuiInterface(const pair<double, double> minPos, const pair<double, double> maxPos, const string& filePath1, const string& filePath2)
	: MicroSimulation<T, U>(minPos, maxPos, filePath1, filePath2)
{}

template<class T, class U>
bool GuiInterface<T, U>::run(int const intervals)
{
	/*
	while (*intervals == -1)
	{
		sleep(1);
	}
	std::cout << "Start" << std::endl;
	std::pair<int,int> pos;
	std::pair<double,double> castPos;
	pos = mGui->getBaseCoordinates(true);
	castPos.first = pos.first;
	castPos.second = pos.second;
	MicroSimulation<T, U>::setPlayer1Pos(castPos);

	pos = mGui->getBaseCoordinates(false);
	castPos.first = pos.first;
	castPos.second = pos.second;
	MicroSimulation<T, U>::setPlayer2Pos(castPos);
	*/

	//std::cout << "Positions set" << std::endl;

	/*
	while (MicroSimulation<T,U>::run(intervals) == false)
	{
		if (*intervals == 0)
		{
			sleep(1);
			continue;
		}
	*/
	//	std::cout << "I'm running" << std::endl;
	bool ret = MicroSimulation<T,U>::run(intervals);
	mGui->setUnits(getPlayer1(), getPlayer2());
	return ret;
	//std::cout << "Finished" << std::endl;
}

template<class T, class U>
void GuiInterface<T,U>::run()
{
    MicroSimulation<T,U>::run(false);
}

template<class T, class U>
void GuiInterface<T,U>::clearBothPlayers()
{
	MicroSimulation<T,U>::clearBothPlayers();
}

template<class T, class U>
void GuiInterface<T,U>::reset()
{
    MicroSimulation<T,U>::resetBothPlayers();
}

template<class T, class U>
void GuiInterface<T, U>::setGuiWindow(GuiWindow *gui)
{
	this->mGui = gui;
}

template<class T, class U>
void GuiInterface<T, U>::initPotentialFields()
{
	MicroSimulation<T,U>::initPotentialFields();
}

template<class T, class U>
void GuiInterface<T, U>::initBothPlayers(const vector<string>& a, const vector<string>& b)
{
	MicroSimulation<T,U>::initBothPlayers(a,b);
	MicroSimulation<T,U>::initPotentialFields();
}

template<class T, class U>
double GuiInterface<T, U>::compareCosts()
{
	return MicroSimulation<T,U>::compareCosts();
}

template<class T, class U>
PlayerState<T>const& GuiInterface<T,U>::getPlayer1() const
{
	return MicroSimulation<T,U>::getPlayer1();
}

template<class T, class U>
PlayerState<U>const& GuiInterface<T,U>::getPlayer2() const
{
	return MicroSimulation<T,U>::getPlayer2();
}

template<class T, class U>
SimulationResult GuiInterface<T, U>::getPlayer1Result()
{
	return MicroSimulation<T,U>::getPlayer1().calculateResult();
}

template<class T, class U>
SimulationResult GuiInterface<T, U>::getPlayer2Result()
{
	return MicroSimulation<T,U>::getPlayer2().calculateResult();
}

template<class T, class U>
void GuiInterface<T,U>::setPositions(std::pair<int,int> &pl1, std::pair<int,int> &pl2)
{
	std::pair<double,double> castPos;
	castPos.first = pl1.first;
	castPos.second = pl1.second;
	MicroSimulation<T, U>::setPlayer1Pos(castPos);

	castPos.first = pl2.first;
	castPos.second = pl2.second;
	MicroSimulation<T, U>::setPlayer2Pos(castPos);
}

template<class T, class U>
int GuiInterface<T, U>::compareStats()
{
	return MicroSimulation<T,U>::compareStats();
}

template<class T, class U>
void GuiInterface<T, U>::collectGarbage()
{
	MicroSimulation<T,U>::collectGarbage();
}

template<class T, class U>
void GuiInterface<T, U>::setTracking(bool const tracking)
{
    MicroSimulation<T,U>::setTracking(tracking);
}

template<class T, class U>
void GuiInterface<T, U>::setTracking(bool const tracking, size_t const steps)
{
    MicroSimulation<T,U>::setTracking(tracking, steps);
}
template<class T, class U>
void GuiInterface<T, U>::setTimeSteps(size_t steps)
{
    MicroSimulation<T,U>::setTimeSteps(steps);
}

template<class T, class U>
void GuiInterface<T, U>::clearUnitPaths()
{
    MicroSimulation<T,U>::clearUnitPaths();
}

