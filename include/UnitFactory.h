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

    //TODO save type in actual object

    void create(const string& name, PlayerState<Race>& pl)
	{
        if (name == "ZergUnit" || name == "TerranUnit" || name == "ProtossUnit")
		{
			return;
		}
		if (name == Race::nameList[0])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT0 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(0);

            pl.unitList0.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList0.back().first);
		}
		else if (name == Race::nameList[1])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT1 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(1);

            pl.unitList1.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList1.back().first);
		}
		else if (name == Race::nameList[2])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT2 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(2);

            pl.unitList2.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList2.back().first);


		}
		else if (name == Race::nameList[3])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT3 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(3);

            pl.unitList3.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList3.back().first);
		}
		else if (name == Race::nameList[4])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT4 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(4);

            pl.unitList4.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList4.back().first);
        }
		else if (name == Race::nameList[5])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT5 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(5);

            pl.unitList5.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList5.back().first);
		}
		else if (name == Race::nameList[6])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT6 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(6);

            pl.unitList6.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList6.back().first);
		}
		else if (name == Race::nameList[7])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT7 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(7);

            pl.unitList7.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList7.back().first);
		}
		else if (name == Race::nameList[8])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT8 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(8);

            pl.unitList8.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList8.back().first);
		}
		else if (name == Race::nameList[9])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT9 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(9);

            pl.unitList9.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList9.back().first);
		}
		else if (name == Race::nameList[10])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT10 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(10);

            pl.unitList10.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList10.back().first);
		}
		else if (name == Race::nameList[11])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT11 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(11);

            pl.unitList11.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList11.back().first);
		}
		else if (name == Race::nameList[12])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT12 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(12);

            pl.unitList12.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList12.back().first);
		}
		else if (name == Race::nameList[13])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT13 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(13);

            pl.unitList13.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList13.back().first);
		}
		else if (name == Race::nameList[14])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT14 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(14);

            pl.unitList14.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList14.back().first);
		}
		else if (name == Race::nameList[15])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT15 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(15);

            pl.unitList15.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList15.back().first);
		}
		else if (name == Race::nameList[16])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT16 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(16);

            pl.unitList16.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList16.back().first);
		}
		else if (name == Race::nameList[17])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            typename Race::UT17 unit;
            unit.setName(name);
            unit.setStats(mUmap[name]);
            unit.setIdentifier(17);

            pl.unitList17.emplace_back(unit,it);
            pl.unitList.back() = &(pl.unitList17.back().first);
		}
	}
    void create(const vector<string>& names, PlayerState<Race>& pl)
	{
		for (const string& name : names)
		{
            create(name, pl);
		}
	}
};

#endif
