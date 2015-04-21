#ifndef _UNITFACTORY_H_
#define _UNITFACTORY_H_

#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <utility>
#include <typeinfo>

#include "Race.h"
#include "Unit.h"
#include "PlayerState.h"


using std::string;
using std::vector;
using std::list;
using std::unordered_map;
using std::pair;
using std::type_info;

template <class Race> 
class UnitFactory : public Race 
{
private:
	unordered_map<string, UnitStats> mUmap;
	
public:
	UnitFactory()
	{}

	UnitFactory(const unordered_map<string, UnitStats> &umap)
        : mUmap(umap)
	{}

    UnitFactory(const UnitFactory<Race>& unitFactory)
        : mUmap(unitFactory.getUmap())
    {}
    
    unordered_map<string,UnitStats> const& getUmap() const
    {
        return mUmap;
    }

	void setUmap(const unordered_map<string, UnitStats> &umap)
	{
		mUmap = umap;
	}

	bool isUmapEmpty() const
	{
		return mUmap.empty();
	}

    template<typename T>
    void createUnit(const string& name, vector<T>& unitListX, int x)
    {
        T unit;
        unit.setName(name);
        unit.setStats(mUmap[name]);
        unit.setIdentifier(x);
        unitListX.push_back (std::move(unit));
    }

    template<typename T>
    void setPointer(PlayerState<Race>& pl, vector<T>& unitListX)
    {
        pl.unitList.reserve (pl.unitList.size() + unitListX.size());
        for(auto& unit : unitListX)
        {
            pl.unitList.emplace_back();
            pl.unitList.back() = &unit;
        }
    }

    void create(const string& name, PlayerState<Race>& pl)
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
    void create(const vector<string>& names, PlayerState<Race>& pl)
	{
		for (const string& name : names)
		{
            create(name, pl);
		}
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
