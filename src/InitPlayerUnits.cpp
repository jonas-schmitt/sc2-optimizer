#include<unordered_map>
#include<string>
#include<algorithm>
#include<fstream>
#include "../include/InitPlayerUnits.h"


template<class Race>
InitPlayerUnits<Race>::InitPlayerUnits()
{}

template<class Race>
InitPlayerUnits<Race>::InitPlayerUnits(const std::string& path) : mFactory(path)
{
    if (path.empty())
	{
        throw std::invalid_argument("InitPlayerUnits::InitPlayerUnits(const std::string&): The file path passed as argument is an empty string.");
	}
    mPath = path;
	readStats();
}
template<class Race>
std::vector<std::string>& InitPlayerUnits<Race>::split(const std::string &s, char delim, std::vector<std::string> &tokens) {
    std::stringstream ss(s);
    std::string token;
    while (getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

template<class Race>
std::vector<std::string> InitPlayerUnits<Race>::split(std::string const &s, char delim) {
    std::vector<std::string> tokens;
    split(s, delim, tokens);
    return tokens;
}

template<class Race>
void InitPlayerUnits<Race>::setPath(const std::string& path)
{
    if (path.empty())
	{
        throw std::invalid_argument("InitPlayerUnits::setpath(const std::string&): The file path passed as argument is an empty string.");
	}
    mPath = path;
    mFactory.setPath (path);
}

template <class Race>
std::string InitPlayerUnits<Race>::getPath() const
{
    return mPath;
}

// Obtain the stats for all units of the player's race
template <class Race>
void InitPlayerUnits<Race>::readStats()
{

    std::ifstream file;
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
        std::string str;
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
        std::vector<std::string> attributes = split(str, ',');
        for(std::string const& attr : attributes)
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
        std::sort(stats.attributes.begin (), stats.attributes.end ());
        if(!(stream
             >> stats.groundAttack
             >> stats.gaUpgrade
             >> stats.airAttack
             >> stats.aaUpgrade))
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
            std::vector<std::string> bonusStrs = split(str, ',');
            for(std::string const& bonusStr : bonusStrs)
            {
                Bonus bonus;
                std::vector<std::string> tokens = split(bonusStr, ':');

                std::stringstream ss(tokens.back ());
                if(!(ss >> bonus.base))
                {
                    bonus.base = 0.f;
                }

                std::vector<std::string> attributes = split(tokens.front(), '&');
                for(std::string const& attr : attributes)
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
                std::sort(bonus.attributes.begin (), bonus.attributes.end());
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
            std::vector<std::string> bonusStrs = split(str,',');
            for(std::string const& bonusStr : bonusStrs)
            {
                Bonus bonus;
                std::vector<std::string> tokens = split(bonusStr, ':');
                std::vector<std::string> attributes = split(tokens.front(), '&');
                for(std::string const& attr : attributes)
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
                std::sort(bonus.attributes.begin (), bonus.attributes.end());
                bool matched = false;
                for(Bonus& b : stats.bonuses)
                {
                    // if there is already a bonus present, just add the upgrade
                    if(bonus.attributes.size () == b.attributes.size ()
                            && std::mismatch(bonus.attributes.begin(), bonus.attributes.end(), b.attributes.begin()).first == bonus.attributes.end ())
                    {
                        std::stringstream ss(*tokens.rbegin ());
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
                    std::stringstream ss(*tokens.rbegin ());
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
void InitPlayerUnits<Race>::init(const std::vector<std::string> &unitVec, const std::string& path, PlayerState<Race>& pl)
{
    if (path.empty())
	{
        throw std::invalid_argument("InitPlayerUnits::init(const std::vector<std::string>&, const std::string&): The file path passed as argument is an empty string.");
	}
	if (unitVec.empty())
	{
        throw std::invalid_argument("InitPlayerUnits::init(const std::vector<std::string>&, const std::string&): The vector containing the names of all units to create is empty");
	}
    // Set the path to the directory containing the unit stats
    setPath(path);
    // Obtain the statistics
	readStats();

    // Create all units
    mFactory.create(unitVec, pl);

    for(auto unit : pl.unitList)
    {
        unit->setPosLimits(pl.minPos, pl.maxPos);
    }
    pl.unitCount = pl.unitList.size ();

    size_t pos = 0;

    // Create the combined chromosome of all units of all players and set the chromosome for each unit accordingly
    if(!pl.unitList0.empty())
    {
        pl.NGenes += pl.unitList0.front().mNGenes;
        for(auto & unit : pl.unitList0)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList1.empty())
    {
        pl.NGenes += pl.unitList1.front().mNGenes;
        for(auto & unit : pl.unitList1)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }

    if(!pl.unitList2.empty())
    {
        pl.NGenes += pl.unitList2.front().mNGenes;
        for(auto & unit : pl.unitList2)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }

    if(!pl.unitList3.empty())
    {
        pl.NGenes += pl.unitList3.front().mNGenes;
        for(auto & unit : pl.unitList3)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }

    if(!pl.unitList4.empty())
    {
        pl.NGenes += pl.unitList4.front().mNGenes;
        for(auto & unit : pl.unitList4)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }

    if(!pl.unitList5.empty())
    {
        pl.NGenes += pl.unitList5.front().mNGenes;
        for(auto & unit : pl.unitList5)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }

    if(!pl.unitList6.empty())
    {
        pl.NGenes += pl.unitList6.front().mNGenes;
        for(auto & unit : pl.unitList6)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }

    if(!pl.unitList7.empty())
    {
        pl.NGenes += pl.unitList7.front().mNGenes;
        for(auto & unit : pl.unitList7)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList8.empty())
    {
        pl.NGenes += pl.unitList8.front().mNGenes;
        for(auto & unit : pl.unitList8)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList9.empty())
    {
        pl.NGenes += pl.unitList9.front().mNGenes;
        for(auto & unit : pl.unitList9)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList10.empty())
    {
        pl.NGenes += pl.unitList10.front().mNGenes;
        for(auto & unit : pl.unitList10)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList11.empty())
    {
        pl.NGenes += pl.unitList11.front().mNGenes;
        for(auto & unit : pl.unitList11)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList12.empty())
    {
        pl.NGenes += pl.unitList12.front().mNGenes;
        for(auto & unit : pl.unitList12)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList13.empty())
    {
        pl.NGenes += pl.unitList13.front().mNGenes;
        for(auto & unit : pl.unitList13)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList14.empty())
    {
        pl.NGenes += pl.unitList14.front().mNGenes;
        for(auto & unit : pl.unitList14)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList15.empty())
    {
        pl.NGenes += pl.unitList15.front().mNGenes;
        for(auto & unit : pl.unitList15)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList16.empty())
    {
        pl.NGenes += pl.unitList16.front().mNGenes;
        for(auto & unit : pl.unitList16)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }
    if(!pl.unitList17.empty())
    {
        pl.NGenes += pl.unitList17.front().mNGenes;
        for(auto & unit : pl.unitList17)
        {
            unit.setChromosomeStartPosition(pos);
        }
        pos = pl.NGenes;
    }


    pl.potentialList.emplace_back(pl.minPos,funcMinX);
    pl.potentialList.emplace_back(pl.maxPos,funcMaxX);
    pl.potentialList.emplace_back(pl.minPos,funcMinY);
    pl.potentialList.emplace_back(pl.maxPos,funcMaxY);

    // The following is required for the dynamic adjustment of the playground size and the initial placement of the units
    // Compute the maximal size
    auto cmp_size = [] (typename Race::RUT *lhs, typename Race::RUT *rhs)
    {
        return lhs->getSize() < rhs->getSize();
    };
    auto it = std::max_element(pl.unitList.begin(), pl.unitList.end(), cmp_size);
    pl.maxUnitSize = (*it)->getSize();

    auto cmp_speed = [] (typename Race::RUT *lhs, typename Race::RUT *rhs)
    {
        return lhs->getSpeed() < rhs->getSpeed();
    };
    std::sort(pl.unitList.begin(), pl.unitList.end(), cmp_speed);

    pl.adjustActionsPerUnit();
}


