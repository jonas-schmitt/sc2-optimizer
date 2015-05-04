#include<unordered_map>
#include<string>
#include<algorithm>
#include<fstream>
#include "../include/InitPlayerUnits.h"
using std::string;
using std::stringstream;
using std::vector;
using std::sort;
using std::ifstream;

template<class Race>
InitPlayerUnits<Race>::InitPlayerUnits()
{}

template<class Race>
InitPlayerUnits<Race>::InitPlayerUnits(const std::string& path) : mFactory(path)
{
    if (path.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::InitPlayerUnits(const string&): The file path passed as argument is an empty string.");
	}
    mPath = path;
	readStats();
}
template<class Race>
vector<string>& InitPlayerUnits<Race>::split(const string &s, char delim, vector<string> &tokens) {
    stringstream ss(s);
    string token;
    while (getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

template<class Race>
vector<string> InitPlayerUnits<Race>::split(string const &s, char delim) {
    vector<string> tokens;
    split(s, delim, tokens);
    return tokens;
}

template<class Race>
void InitPlayerUnits<Race>::setPath(const std::string& path)
{
    if (path.empty())
	{
        throw std::invalid_argument("InitPlayerUnits::setpath(const string&): The file path passed as argument is an empty string.");
	}
    mPath = path;
    mFactory.setPath (path);
}

template <class Race>
std::string InitPlayerUnits<Race>::getPath() const
{
    return mPath;
}


template <class Race>
void InitPlayerUnits<Race>::readStats()
{

    ifstream file;
    file.open (mPath+"/stats.txt");
	std::string line;
	std::string name;

    std::unordered_map<std::string, UnitStats> statMap;
    while (std::getline(file, line))
	{
		std::stringstream stream(line);
		if (!(stream >> name) || name == "Name" || name.empty())
		{
			continue;

		}
        UnitStats stats;
        string str;
        if (!(stream
              >> stats.minerals
              >> stats.gas
              >> stats.size
              >> stats.armor
              >> stats.armorUpgradeBonus
              >> stats.health
              >> stats.shield))
        {

            continue;
        }



        if(!(stream >> str))
        {
            continue;
        }

        // parse attributes
        vector<string> attributes = split(str, ',');
        for(string const& attr : attributes)
        {
            if(attr == "L")
            {
                stats.attributes.push_back(Attribute::light);
            }
            else if(attr == "A")
            {
                stats.attributes.push_back (Attribute::armored);
            }
            else if(attr == "B")
            {
                stats.attributes.push_back (Attribute::biological);
            }
            else if(attr == "M")
            {
                stats.attributes.push_back(Attribute::mechanical);
            }
            else if(attr == "P")
            {
                stats.attributes.push_back (Attribute::psyonic);
            }
            else if(attr == "Ma")
            {
                stats.attributes.push_back (Attribute::massive);
            }
            else if(attr == "Air")
            {
                stats.attributes.push_back (Attribute::air);
                stats.airUnit = true;
            }
        }
        sort(stats.attributes.begin (), stats.attributes.end ());
        if(!(stream
             >> stats.groundAttack
             >> stats.gaUpgradeBonus
             >> stats.airAttack
             >> stats.aaUpgradeBonus))
        {
            continue;
        }

        if(!(stream >> str))
        {
            continue;
        }

        //parse bonuses
        if(str != "-")
        {
            vector<string> bonusStrs = split(str, ',');
            for(string const& bonusStr : bonusStrs)
            {
                Bonus bonus;
                vector<string> tokens = split(bonusStr, ':');

                stringstream ss(tokens.back ());
                if(!(ss >> bonus.base))
                {
                    bonus.base = 0.f;
                }

                vector<string> attributes = split(tokens.front(), '&');
                for(string const& attr : attributes)
                {

                    if(attr == "L")
                    {
                        bonus.attributes.push_back(Attribute::light);
                    }
                    else if(attr == "A")
                    {
                        bonus.attributes.push_back (Attribute::armored);
                    }
                    else if(attr == "B")
                    {
                        bonus.attributes.push_back (Attribute::biological);
                    }
                    else if(attr == "M")
                    {
                        bonus.attributes.push_back(Attribute::mechanical);
                    }
                    else if(attr == "P")
                    {
                        bonus.attributes.push_back (Attribute::psyonic);
                    }
                    else if(attr == "Ma")
                    {
                        bonus.attributes.push_back (Attribute::massive);
                    }
                    else if(attr == "Air")
                    {
                        bonus.attributes.push_back (Attribute::air);
                    }
                }
                sort(bonus.attributes.begin (), bonus.attributes.end());
                stats.bonuses.push_back(bonus);
            }
        }
        //parse bonus upgrades
        if(!(stream >> str))
        {
            continue;
        }
        if(str != "-")
        {
            vector<string> bonusStrs = split(str,',');
            for(string const& bonusStr : bonusStrs)
            {
                Bonus bonus;
                vector<string> tokens = split(bonusStr, ':');
                vector<string> attributes = split(tokens.front(), '&');
                for(string const& attr : attributes)
                {
                    if(attr == "L")
                    {
                        bonus.attributes.push_back(Attribute::light);
                    }
                    else if(attr == "A")
                    {
                        bonus.attributes.push_back (Attribute::armored);
                    }
                    else if(attr == "B")
                    {
                        bonus.attributes.push_back (Attribute::biological);
                    }
                    else if(attr == "M")
                    {
                        bonus.attributes.push_back(Attribute::mechanical);
                    }
                    else if(attr == "P")
                    {
                        bonus.attributes.push_back (Attribute::psyonic);
                    }
                    else if(attr == "Ma")
                    {
                        bonus.attributes.push_back (Attribute::massive);
                    }
                    else if(attr == "Air")
                    {
                        bonus.attributes.push_back (Attribute::air);
                    }
                }
                sort(bonus.attributes.begin (), bonus.attributes.end());
                bool matched = false;
                for(Bonus& b : stats.bonuses)
                {
                    // if there is already a bonus present, just add the upgrade
                    if(bonus.attributes.size () == b.attributes.size ()
                            && std::mismatch(bonus.attributes.begin(), bonus.attributes.end(), b.attributes.begin()).first == bonus.attributes.end ())
                    {
                        stringstream ss(*tokens.rbegin ());
                        if(!(ss >> b.upgrade))
                        {
                            b.upgrade = 0.f;
                        }
                        matched = true;
                        break;
                    }
                }
                // if there is no base but just an upgrade, add this as an additional bonus
                if(!matched)
                {
                    stringstream ss(*tokens.rbegin ());
                    if(!(ss >> bonus.upgrade))
                    {
                        bonus.upgrade = 0.f;
                    }
                    stats.bonuses.push_back (bonus);
                }
            }
        }
        double tmp;

        // the cooldown needs to be rounded to milliseconds
        if(!(stream >> tmp))
        {
            continue;
        }
        stats.gaCooldown = static_cast<int>(std::round(1000*tmp));
        if(!(stream >> tmp))
        {
            continue;
        }
        stats.aaCooldown = static_cast<int>(std::round(1000*tmp));
        if(!(stream >> stats.speed))
        {
            continue;
        }
        if(!(stream
             >> stats.groundRange
             >> stats.airRange
             >> stats.energy
             >> stats.maxEnergy
             >> stats.creepMultiplier))
        {
            continue;
        }
        stats.maxHealth = stats.health;
        stats.maxShield = stats.shield;
        stats.sumMaxHealthAndShield = stats.maxHealth + stats.maxShield;
        statMap[name] = stats;

	}
    mFactory.setHashMap(statMap);
    file.close ();
}


template <class Race>
void InitPlayerUnits<Race>::init(const std::vector<std::string> &unitVec, PlayerState<Race>& pl)
{
    if (mPath.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::init(const vector<string>&): The file path currently stored is an empty string.");
	}
	if (unitVec.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::init(const vector<string>&): The vector containing the names of all units to create is empty");
	}
    if (mFactory.isHashMapEmpty())
	{
		readStats();
	}
    mFactory.create(unitVec, pl);
    for(auto unit : pl.unitList)
    {
        unit->setMinPos(pl.minPos);
        unit->setMaxPos(pl.maxPos);
    }
    pl.unitCount = pl.unitList.size ();

    pl.potentialList.emplace_back(pl.minPos,funcMinX);
    pl.potentialList.emplace_back(pl.maxPos,funcMaxX);
    pl.potentialList.emplace_back(pl.minPos,funcMinY);
    pl.potentialList.emplace_back(pl.maxPos,funcMaxY);
}

template <class Race>
void InitPlayerUnits<Race>::init(const std::vector<std::string> &unitVec, const std::string& path, PlayerState<Race>& pl)
{
    if (path.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::init(const vector<string>&, const string&): The file path passed as argument is an empty string.");
	}
	if (unitVec.empty())
	{
		throw std::invalid_argument("InitPlayerUnits::init(const vector<string>&, const string&): The vector containing the names of all units to create is empty");
	}
    setPath(path);
	readStats();
    mFactory.create(unitVec, pl);

    for(auto unit : pl.unitList)
    {
        unit->setMinPos(pl.minPos);
        unit->setMaxPos(pl.maxPos);
    }
    pl.unitCount = pl.unitList.size ();

    pl.potentialList.emplace_back(pl.minPos,funcMinX);
    pl.potentialList.emplace_back(pl.maxPos,funcMaxX);
    pl.potentialList.emplace_back(pl.minPos,funcMinY);
    pl.potentialList.emplace_back(pl.maxPos,funcMaxY);
}


