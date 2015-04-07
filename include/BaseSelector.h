#ifndef _BASE_SELECTOR_H
#define _BASE_SELECTOR_H

#include <iostream>
#include <unistd.h>

#include "MicroSimulation.h"
#include "Unit.h"

class BaseSelector
{
	public:
        virtual bool run(size_t intervals)=0;
        virtual void run()=0;
        virtual void reset()=0;

		virtual void initBothPlayers(const vector<string>&, const vector<string>&)=0;
		virtual void initPotentialFields()=0;
		virtual double compareCosts()=0;

		virtual SimulationResult getPlayer1Result()=0;
		virtual SimulationResult getPlayer2Result()=0;

		virtual void setTracking(const bool, const size_t)=0;
		virtual void setTracking(const bool)=0;

        virtual void setTimeSteps(size_t)=0;
		
		virtual void setPlayer1Genes(UnitGenes const& genes)=0;
		virtual void setPlayer2Genes(UnitGenes const& genes)=0;

		virtual void getPotentialField()=0;

		virtual void clearBothPlayers()=0;
		virtual void setPositions(std::pair<int, int> &, std::pair<int, int> &)=0; 

		virtual int compareStats()=0;
		virtual void collectGarbage()=0;
};

#endif
