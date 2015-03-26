#ifndef _RACE_SELECTOR_H_
#define _RACE_SELECTOR_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>

#include "../include/Race.h"
#include "../include/UnitFactory.h"

bool findRace(std::vector<std::string> &, std::string &);
void minimizeInputFile(std::vector<std::string> &, std::string path);

#endif
