#ifndef _INITPLAYERUNITS_H_
#define _INITPLAYERUNITS_H_

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <sstream>

#include "DataReader.h"
#include "PlayerState.h"
#include "Race.h"
#include "Unit.h"
#include "UnitFactory.h"
#include "Utilities.h"
using std::string;


template<class Race>
class InitPlayerUnits final : public Race
{

private:
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


    UnitFactory<Race> mFactory;
    std::string mPath;
    vector<string>& split(const string &s, char delim, vector<string> &tokens);
    vector<string> split(string const &s, char delim);
public:
	InitPlayerUnits();
    InitPlayerUnits(const std::string& path);
    void readStats();
    std::string getPath() const;
    void setPath(const std::string& path);
    void init(const std::vector<std::string> &unitVec, PlayerState<Race>& pl);
    void init (const std::vector<std::string>& unitVec, const std::string& path, PlayerState<Race>& pl);
	
};

#endif

