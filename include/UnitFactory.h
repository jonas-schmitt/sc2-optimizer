#ifndef _UNITFACTORY_H_
#define _UNITFACTORY_H_

#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Race.h"
#include "Unit.h"
#include "PlayerState.h"

// Factory class for creating units
template <class Race> 
class UnitFactory final : public Race
{
private:
    std::unordered_map<std::string, UnitStats> mHashMap;
    std::string mPath;
	
public:
	UnitFactory()
	{}

    UnitFactory(std::string const& path)
        : mPath(path)
    {}

    UnitFactory(std::unordered_map<std::string, UnitStats> const &hashMap)
        : mHashMap(hashMap)
	{}

    UnitFactory(std::unordered_map<std::string, UnitStats> const &hashMap, std::string const& path)
        : mHashMap(hashMap), mPath(path)
    {}

    UnitFactory(const UnitFactory<Race>& unitFactory)
        : mHashMap(unitFactory.getHashMap()), mPath(mPath)
    {}

    std::string getPath() const
    {
        return mPath;
    }

    void setPath(const std::string &path)
    {
        mPath = path;
    }

    std::unordered_map<std::string,UnitStats> const& getHashMap() const
    {
        return mHashMap;
    }

    void setHashMap(const std::unordered_map<std::string, UnitStats> &hashMap)
	{
        mHashMap = hashMap;
	}


    bool isHashMapEmpty() const
	{
        return mHashMap.empty();
	}

    // Create a list of units of the same type
    template<typename T>
    void createUnit(const std::string& name, std::vector<T>& unitListX, int x)
    {
        T unit;
        unit.setName(name);
        unit.setStats(mHashMap[name]);
        std::ifstream file;
        file.open(mPath+"/upgrades/"+name+".txt");
        std::string buf;
        std::getline(file, buf);
        std::getline(file, buf);
        int flag;
        std::stringstream ss(buf);
        std::vector<int> flags;
        while(ss >> flag)
        {
            flags.push_back (flag);
        }
        file.close ();
        if(flags.size () < 2)
        {
            flags.resize(2);
        }
        // this method sets the upgrades according to the flags that are passed
        unit.initUpgrades(flags);

        unit.setIdentifier(x);
        unitListX.push_back (std::move(unit));
    }

    // Set pointers to a list of units of the same type
    template<typename T>
    void setPointer(PlayerState<Race>& pl, std::vector<T>& unitListX)
    {
        pl.unitList.reserve (pl.unitList.size() + unitListX.size());
        for(auto& unit : unitListX)
        {
            pl.unitList.emplace_back();
            pl.unitList.back() = &unit;
        }
    }

    // Create a unit
    // name: name of the unit
    void create(const std::string& name, PlayerState<Race>& pl)
	{
        if (name == "ZergUnit" || name == "TerranUnit" || name == "ProtossUnit")
		{
            return;
        }
        if (name == Race::nameList[0])
        {
            createUnit(name, pl.unitList0, 0);
        }
        else if (name == Race::nameList[1])
        {
            createUnit(name, pl.unitList1, 1);
        }
        else if (name == Race::nameList[2])
        {
            createUnit(name, pl.unitList2, 2);
        }
        else if (name == Race::nameList[3])
        {
            createUnit(name, pl.unitList3, 3);
        }
        else if (name == Race::nameList[4])
        {
            createUnit(name, pl.unitList4, 4);
        }
        else if (name == Race::nameList[5])
        {
            createUnit(name, pl.unitList5, 5);
        }
        else if (name == Race::nameList[6])
        {
            createUnit(name, pl.unitList6, 6);
        }
        else if (name == Race::nameList[7])
        {
            createUnit(name, pl.unitList7, 7);
        }
        else if (name == Race::nameList[8])
        {
            createUnit(name, pl.unitList8, 8);
        }
        else if (name == Race::nameList[9])
        {
            createUnit(name, pl.unitList9, 9);
        }
        else if (name == Race::nameList[10])
        {
            createUnit(name, pl.unitList11, 10);
        }
        else if (name == Race::nameList[11])
        {
            createUnit(name, pl.unitList11, 11);
        }
        else if (name == Race::nameList[12])
        {
            createUnit(name, pl.unitList12, 12);
        }
        else if (name == Race::nameList[13])
        {
            createUnit(name, pl.unitList13, 13);
        }
        else if (name == Race::nameList[14])
        {
            createUnit(name, pl.unitList14, 14);
        }
        else if (name == Race::nameList[15])
        {
            createUnit(name, pl.unitList15, 15);
        }
        else if (name == Race::nameList[16])
        {
            createUnit(name, pl.unitList16, 16);
        }
        else if (name == Race::nameList[17])
        {
            createUnit(name, pl.unitList17, 17);
		}
	}

    // Create a number of units and store them in the PlayerState
    // names: Names of the units
    // pl: State of the player for whom the units should be created
    void create(const std::vector<std::string>& names, PlayerState<Race>& pl)
	{
        for (const std::string& name : names)
		{
            create(name, pl);
		}

        // Set the pointers in pl accordingly
        setPointer(pl, pl.unitList0);
        setPointer(pl, pl.unitList1);
        setPointer(pl, pl.unitList2);
        setPointer(pl, pl.unitList3);
        setPointer(pl, pl.unitList4);
        setPointer(pl, pl.unitList5);
        setPointer(pl, pl.unitList6);
        setPointer(pl, pl.unitList7);
        setPointer(pl, pl.unitList8);
        setPointer(pl, pl.unitList9);
        setPointer(pl, pl.unitList10);
        setPointer(pl, pl.unitList11);
        setPointer(pl, pl.unitList12);
        setPointer(pl, pl.unitList13);
        setPointer(pl, pl.unitList14);
        setPointer(pl, pl.unitList15);
        setPointer(pl, pl.unitList16);
        setPointer(pl, pl.unitList17);

	}
};

#endif
