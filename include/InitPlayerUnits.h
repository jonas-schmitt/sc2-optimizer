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
class InitPlayerUnits : public Race
{

private:
    std::function<Vec2D(Vec2D const& pos,typename Race::BUT const& unit)> funcMinX = [](Vec2D const& pos,typename Race::BUT const& unit)
    {
        double const dist = fabs(pos.x-unit.getX());
        if(dist < unit.getSpeed())
        {
            return Vec2D(MAX*1000,0);
        }
        double const z = 1/pow(dist,3)*MAX*1000;

        return Vec2D(z,0);
    };
    std::function<Vec2D(Vec2D const& pos,typename Race::BUT const& unit)> funcMinY = [](Vec2D const& pos,typename Race::BUT const& unit)
    {
        double const dist = fabs(pos.y-unit.getY());
        if(dist < unit.getSpeed())
        {
            return Vec2D(0,MAX*1000);
        }

        double const z = 1/pow(dist,3)*MAX*1000;
        return Vec2D(0,z);
    };
    std::function<Vec2D(Vec2D const& pos,typename Race::BUT const& unit)> funcMaxX = [](Vec2D const& pos,typename Race::BUT const& unit)
    {
        double const dist = fabs(pos.x-unit.getX());
        if(dist < unit.getSpeed())
        {
            return Vec2D(-MAX*1000,0);
        }

        double const z = 1/pow(dist,3)*MAX*1000;
        return Vec2D(-z,0);
    };
    std::function<Vec2D(Vec2D const& pos,typename Race::BUT const& unit)> funcMaxY = [](Vec2D const& pos,typename Race::BUT const& unit)
    {
        double const dist = fabs(pos.y-unit.getY());
        if(dist < unit.getSpeed())
        {
            return Vec2D(0,-MAX*1000);
        }

        double const z = 1/pow(dist,3)*MAX*1000;
        return Vec2D(0,-z);
    };


    UnitFactory<Race> mFactory;
	std::string mFilePath;
    vector<string>& split(const string &s, char delim, vector<string> &tokens);
    vector<string> split(string const &s, char delim);
public:
	InitPlayerUnits();
	InitPlayerUnits(const std::string& filePath);
    void readStats();
	std::string getFilePath() const;
	void setFilePath(const std::string& filePath);
    void init(const std::vector<std::string> &unitVec, PlayerState<Race>& pl);
    void init (const std::vector<std::string>& unitVec, const std::string& filePath, PlayerState<Race>& pl);
	
};

#endif

