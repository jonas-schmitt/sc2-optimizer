#include<unordered_map>
#include<string>
#include "../include/InitPlayerUnits.h"
using std::string;

template<class Race>
InitPlayerUnits<Race>::InitPlayerUnits()
{}

template<class Race>
InitPlayerUnits<Race>::InitPlayerUnits(const std::string& filePath)
{
	if (filePath.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::InitPlayerUnits(const string&): The file path passed as argument is an empty string.");
	}
	mFilePath = filePath;
	readStats();
}

template<class Race>
void InitPlayerUnits<Race>::setFilePath(const std::string& filePath)
{
	if (filePath.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::setFilePath(const string&): The file path passed as argument is an empty string.");
	}
	mFilePath = filePath;
}

template <class Race>
std::string InitPlayerUnits<Race>::getFilePath() const
{
	return mFilePath;
}


template <class Race>
void InitPlayerUnits<Race>::readStats()
{

	DataReader reader(mFilePath);
	std::string line;
	UnitStats stats;
	std::string name;

    std::unordered_map<std::string, UnitStats> statMap;

	while (!((line = reader.getLine()).empty()))
	{
		std::stringstream stream(line);
		if (!(stream >> name) || name == "Name" || name.empty())
		{
			continue;
		}
        string attrStr;
		if (!(stream
			>> stats.minerals
			>> stats.gas
            >> stats.size
            >> stats.armor
            >> stats.armorUpgrade
            >> stats.health
            >> stats.shield))
		{
			continue;
		}
        if(!(stream >> attrStr))
        {
            continue;
        }


        stats.maxHealth = stats.health;
        stats.maxShield = stats.shield;
        stats.maxEnergy = stats.energy;
        statMap[name] = stats;
	}
    mFactory.setUmap(statMap);
}


template <class Race>
void InitPlayerUnits<Race>::init(const std::vector<std::string> &unitVec, PlayerState<Race>& pl)
{
	if (mFilePath.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::init(const vector<string>&): The file path currently stored is an empty string.");
	}
	if (unitVec.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::init(const vector<string>&): The vector containing the names of all units to create is empty");
	}
	if (mFactory.isUmapEmpty())
	{
		readStats();
	}
    mFactory.create(unitVec, pl);
    for(auto unit : pl.unitList)
    {
        unit->setMinPos(pl.minPos);
        unit->setMaxPos(pl.maxPos);
    }
}

template <class Race>
void InitPlayerUnits<Race>::init(const std::vector<std::string> &unitVec, const std::string& filePath, PlayerState<Race>& pl)
{
	if (filePath.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::init(const vector<string>&, const string&): The file path passed as argument is an empty string.");
	}
	if (unitVec.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::init(const vector<string>&, const string&): The vector containing the names of all units to create is empty");
	}
	setFilePath(filePath);
	readStats();
    mFactory.create(unitVec, pl);
    for(auto unit : pl.unitList)
    {
        unit->setMinPos(pl.minPos);
        unit->setMaxPos(pl.maxPos);
    }
}


