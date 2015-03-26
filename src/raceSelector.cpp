#include "../include/raceSelector.h"

template <class Race>
bool findString(std::string findMe)
{
	Race race;
	for (int i = 0; i < 18; ++i)
	{
		if (findMe == race.nameList[i])
			return true;
	}
	return false;
}

bool findRace(std::vector<std::string> &in, std::string &ret)
{
	bool found = false;
	size_t schlepp = 0;

	while (found == false && schlepp < in.size())
	{
		if (findString<Terran>(in.at(schlepp)) == true)
		{
			ret = "Terran";
		}
		else if (findString<Zerg>(in.at(schlepp)) == true)
		{
			ret = "Zerg";
		}
		else if (findString<Protoss>(in.at(schlepp)) == true)
		{
			ret = "Protoss";
		} else
		{
			++schlepp;
			continue;
		}
		found=true;
		break;
	}
	return found;
}

void minimizeInputFile(std::vector<std::string> &units, std::string path)
{
	std::ifstream in(path.c_str(), std::ifstream::in);
	char tmp[256];
	in.getline(tmp, 256);
	while (in.good())
	{
		std::stringstream ss;
		in.getline(tmp, 256);
		ss << tmp;
		while (ss.good())
		{
			ss >> tmp;
		}
        std::string unitName(tmp);
        if(unitName != "Probe" && unitName != "SCV" && unitName != "Drone" && unitName != "Overlord")
        {
            units.push_back(unitName);
            if(unitName == "Zergling")
            {
                units.push_back(unitName);
            }
        }
	}
}

