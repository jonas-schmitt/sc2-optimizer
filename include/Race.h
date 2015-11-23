#ifndef _RACE_H_
#define _RACE_H_

#include<array>
#include<string>

#include "Unit.h"

// Type definitions for the race types

struct Terran
{
    typedef BaseUnit BUT;
    typedef TerranUnit RUT;
	typedef SCV UT0;
	typedef Marine UT1;
	typedef Marauder UT2;
	typedef Reaper UT3;
	typedef Ghost UT4;
	typedef Hellion UT5;
	typedef Hellbat UT6;
	typedef SiegeTank UT7;
	typedef WidowMine UT8;
	typedef Thor UT9;
	typedef Viking UT10;
	typedef Medivac UT11;
	typedef Raven UT12;
	typedef Banshee UT13;
	typedef Battlecruiser UT14;
    typedef TerranUnit UT15;
    typedef TerranUnit UT16;
    typedef TerranUnit UT17;
    std::array<std::string, 18> nameList{{ "SCV", "Marine", "Marauder", "Reaper", "Ghost", "Hellion", "Hellbat", "SiegeTank", "WidowMine", "Thor",
        "Viking", "Medivac", "Raven", "Banshee", "Battlecruiser", "TerranUnit", "TerranUnit", "TerranUnit" } };

};

struct Zerg
{
    typedef BaseUnit BUT;
    typedef ZergUnit RUT;
	typedef Drone UT0;
	typedef Queen UT1;
	typedef Zergling UT2;
	typedef Baneling UT3;
	typedef Roach UT4;
	typedef Hydralisk UT5;
	typedef Infestor UT6;
	typedef SwarmHost UT7;
	typedef Ultralisk UT8;
	typedef Overseer UT9;
	typedef Mutalisk UT10;
	typedef Corruptor UT11;
	typedef BroodLord UT12;
	typedef Viper UT13;
    typedef ZergUnit UT14;
    typedef ZergUnit UT15;
    typedef ZergUnit UT16;
    typedef ZergUnit UT17;
    std::array<std::string, 18> nameList{{ "Drone", "Queen", "Zergling", "Baneling", "Roach", "Hydralisk", "Infestor", "SwarmHost", "Ultralisk",
        "Overseer", "Mutalisk", "Corruptor", "BroodLord", "Viper", "ZergUnit", "ZergUnit", "ZergUnit", "ZergUnit" } };
    
};

struct Protoss
{
    typedef BaseUnit BUT;
    typedef ProtossUnit RUT;
	typedef Probe UT0;
	typedef Zealot UT1;
	typedef Stalker UT2;
	typedef Sentry UT3;
	typedef HighTemplar UT4;
	typedef DarkTemplar UT5;
	typedef Immortal UT6;
	typedef Colossus UT7;
	typedef Archon UT8;
	typedef Observer UT9;
	typedef WarpPrism UT10;
	typedef Phoenix UT11;
	typedef VoidRay UT12;
	typedef Oracle UT13;
	typedef Carrier UT14;
	typedef Tempest UT15;
	typedef MothershipCore UT16;
	typedef Mothership UT17;
    std::array<std::string, 18> nameList{{ "Probe", "Zealot", "Stalker", "Sentry", "HighTemplar", "DarkTemplar", "Immortal", "Colossus", "Archon",
		"Observer", "WarpPrism", "Phoenix", "VoidRay", "Oracle", "Carrier", "Tempest", "MothershipCore", "Mothership" } };
    
};


#endif
