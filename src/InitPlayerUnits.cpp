#include<unordered_map>
#include<string>
#include<algorithm>
#include "../include/InitPlayerUnits.h"
using std::string;
using std::stringstream;
using std::vector;
using std::sort;

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
        string str;
		if (!(stream
			>> stats.minerals
			>> stats.gas
            >> stats.size
            >> stats.armor
            >> stats.armorUpgrade
            >> stats.health
            >> stats.shield))
		{
            std::cout << "Skip " << name << " " << 1 << std::endl;
			continue;
		}


        if(!(stream >> str))
        {
                        std::cout << "Skip " << name << " " << 2 << std::endl;
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
             >> stats.gaUpgrade
             >> stats.airAttack
             >> stats.aaUpgrade))
        {
                        std::cout << "Skip " << name << " " << 3 << std::endl;
            continue;
        }

        if(!(stream >> str))
        {
                        std::cout << "Skip " << name << " " << 4 << std::endl;
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
                        std::cout << "Skip " << name << " " << 5 << std::endl;
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
        if(!(stream >> stats.gaCooldown))
        {
            std::cout << "Skip " << name << ": gaCooldown" << std::endl;
            continue;
        }
        if(!(stream >> stats.aaCooldown))
        {
            std::cout << "Skip " << name << ": aaCooldown" << std::endl;
            continue;


        }
        if(!(stream >> stats.speed))
        {
            std::cout << "Skip " << name << ": speed" << std::endl;
            continue;
        }
        if(!(stream
             >> stats.groundRange
             >> stats.airRange))
        {
            std::cout << "Skip " << name << ": range" << std::endl;
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


