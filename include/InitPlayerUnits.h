#ifndef _INITPLAYERUNITS_H_
#define _INITPLAYERUNITS_H_

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <sstream>

#include "DataReader.h"
#include "PlayerState.h"
#include "Race.h"
#include "Unit.h"
#include "UnitFactory.h"
using std::string;


template<class Race>
class InitPlayerUnits : public Race
{

private:
    UnitFactory<Race> mFactory;
	std::string mFilePath;
    vector<string>& split(const string &s, char delim, vector<string> &tokens);
    vector<string> split(string const &s, char delim);
public:
	InitPlayerUnits();
	InitPlayerUnits(const std::string& filePath);
    void readStats();
	std::string getFilePath() const;
	void setFilePath(const std::string& filePath);
    void init(const std::vector<std::string> &unitVec, PlayerState<Race>& pl);
    void init (const std::vector<std::string>& unitVec, const std::string& filePath, PlayerState<Race>& pl);
	
};

#endif

