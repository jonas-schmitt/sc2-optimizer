#ifndef _INITPLAYERUNITS_H_
#define _INITPLAYERUNITS_H_

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <sstream>

#include "PlayerState.h"
#include "Race.h"
#include "Unit.h"
#include "UnitFactory.h"
#include "Utilities.h"


// Class for initialising the units of a certain player

template<class Race>
class InitPlayerUnits final : public Race
{

private:

    // Potential fields that force units to avoid the borders of the playground
    std::function<Vec2D(Vec2D const& pos,typename Race::BUT const& unit)> funcMinX = [](Vec2D const& pos,typename Race::BUT const& unit)
    {
        double const dist = fabs(pos.x-unit.getX());
        if(dist < unit.getMovementUpdateDist() + unit.getSize())
        {
            return Vec2D(1e6,0);
        }
        double const z = pow(dist,-3)*1e6;

        return Vec2D(z,0);
    };
    std::function<Vec2D(Vec2D const& pos,typename Race::BUT const& unit)> funcMinY = [](Vec2D const& pos,typename Race::BUT const& unit)
    {
        double const dist = fabs(pos.y-unit.getY());
        if(dist < unit.getMovementUpdateDist() + unit.getSize())
        {
            return Vec2D(0,1e6);
        }

        double const z = pow(dist,-3)*1e6;
        return Vec2D(0,z);
    };
    std::function<Vec2D(Vec2D const& pos,typename Race::BUT const& unit)> funcMaxX = [](Vec2D const& pos,typename Race::BUT const& unit)
    {
        double const dist = fabs(pos.x-unit.getX());
        if(dist < unit.getMovementUpdateDist() + unit.getSize())
        {
            return Vec2D(-1e6,0);
        }

        double const z = pow(dist,-3)*1e6;
        return Vec2D(-z,0);
    };
    std::function<Vec2D(Vec2D const& pos,typename Race::BUT const& unit)> funcMaxY = [](Vec2D const& pos,typename Race::BUT const& unit)
    {
        double const dist = fabs(pos.y-unit.getY());
        if(dist < unit.getMovementUpdateDist() + unit.getSize())
        {
            return Vec2D(0,-1e6);
        }

        double const z = pow(dist,-3)*1e6;
        return Vec2D(0,-z);
    };


    // Factory for the creation of units
    UnitFactory<Race> mFactory;

    // Path to the files containing the different unit stats
    std::string mPath;

    // Helper functions for parsing
    std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &tokens);
    std::vector<std::string> split(std::string const &s, char delim);

public:
	InitPlayerUnits();
    InitPlayerUnits(const std::string& path);

    // Obtain all statistics for all units of the player's race
    void readStats();
    std::string getPath() const;
    void setPath(const std::string& path);

    // Initialize the state of the player according to a certain build order
    // unitVec: Names of the units to create
    // path: Path to the directory containing the unit stats
    // pl: State of the player
    void init (const std::vector<std::string>& unitVec, const std::string& path, PlayerState<Race>& pl);
	
};

#endif

