#ifndef _UNITFACTORY_H_
#define _UNITFACTORY_H_

#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <utility>

#include "Race.h"
#include "Unit.h"
#include "PlayerState.h"

using std::string;
using std::vector;
using std::list;
using std::unordered_map;
using std::pair;

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
            pl.unitList0.emplace_back(typename Race::UT0(),it);
            pl.unitList0.back().first.setName(Race::nameList[0]);
            pl.unitList0.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList0.back().first);
		}
		else if (name == Race::nameList[1])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList1.emplace_back(typename Race::UT1(),it);
            pl.unitList1.back().first.setName(Race::nameList[1]);
            pl.unitList1.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList1.back().first);
		}
		else if (name == Race::nameList[2])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList2.emplace_back(typename Race::UT2(),it);
            pl.unitList2.back().first.setName(Race::nameList[2]);
            pl.unitList2.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList2.back().first);

		}
		else if (name == Race::nameList[3])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList3.emplace_back(typename Race::UT3(),it);
            pl.unitList3.back().first.setName(Race::nameList[3]);
            pl.unitList3.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList3.back().first);
		}
		else if (name == Race::nameList[4])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList4.emplace_back(typename Race::UT4(),it);
            pl.unitList4.back().first.setName(Race::nameList[4]);
            pl.unitList4.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList4.back().first);
        }
		else if (name == Race::nameList[5])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList5.emplace_back(typename Race::UT5(),it);
            pl.unitList5.back().first.setName(Race::nameList[5]);
            pl.unitList5.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList5.back().first);
		}
		else if (name == Race::nameList[6])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList6.emplace_back(typename Race::UT6(),it);
            pl.unitList6.back().first.setName(Race::nameList[6]);
            pl.unitList6.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList6.back().first);
		}
		else if (name == Race::nameList[7])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList7.emplace_back(typename Race::UT7(),it);
            pl.unitList7.back().first.setName(Race::nameList[7]);
            pl.unitList7.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList7.back().first);
		}
		else if (name == Race::nameList[8])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList8.emplace_back(typename Race::UT8(),it);
            pl.unitList8.back().first.setName(Race::nameList[8]);
            pl.unitList8.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList8.back().first);
		}
		else if (name == Race::nameList[9])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList9.emplace_back(typename Race::UT9(),it);
            pl.unitList9.back().first.setName(Race::nameList[9]);
            pl.unitList9.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList9.back().first);
		}
		else if (name == Race::nameList[10])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList10.emplace_back(typename Race::UT10(),it);
            pl.unitList10.back().first.setName(Race::nameList[10]);
            pl.unitList10.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList10.back().first);
		}
		else if (name == Race::nameList[11])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList11.emplace_back(typename Race::UT11(),it);
            pl.unitList11.back().first.setName(Race::nameList[11]);
            pl.unitList11.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList11.back().first);
		}
		else if (name == Race::nameList[12])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList12.emplace_back(typename Race::UT12(),it);
            pl.unitList12.back().first.setName(Race::nameList[12]);
            pl.unitList12.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList12.back().first);
		}
		else if (name == Race::nameList[13])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList13.emplace_back(typename Race::UT13(),it);
            pl.unitList13.back().first.setName(Race::nameList[13]);
            pl.unitList13.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList13.back().first);
		}
		else if (name == Race::nameList[14])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList14.emplace_back(typename Race::UT14(),it);
            pl.unitList14.back().first.setName(Race::nameList[14]);
            pl.unitList14.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList14.back().first);
		}
		else if (name == Race::nameList[15])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList15.emplace_back(typename Race::UT15(),it);
            pl.unitList15.back().first.setName(Race::nameList[15]);
            pl.unitList15.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList15.back().first);
		}
		else if (name == Race::nameList[16])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList16.emplace_back(typename Race::UT16(),it);
            pl.unitList16.back().first.setName(Race::nameList[16]);
            pl.unitList16.back().first.setStats(mUmap[name]);
            pl.unitList.back() = &(pl.unitList16.back().first);
		}
		else if (name == Race::nameList[17])
		{
            pl.unitList.emplace_back();
            auto it = pl.unitList.end();
            --it;
            pl.unitList17.emplace_back(typename Race::UT17(),it);
            pl.unitList17.back().first.setName(Race::nameList[17]);
            pl.unitList17.back().first.setStats(mUmap[name]);
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
