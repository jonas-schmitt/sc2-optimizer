#ifndef ZERG_H
#define ZERG_H

#include<utility>
#include<functional>
#include<cmath>
#include<list>
#include<chrono>
#include<random>
#include<string>
#include<vector>
#include<typeinfo>
#include<array>
#include<deque>


#include "PlayerState.h"
#include "Utilities.h"
#include "UnitGenes.h"
#include "BaseUnit.h"

using std::pair;
using std::tuple;
using std::function;
using std::pow;
using std::sqrt;
using std::list;
using std::string;
using std::vector;
using std::array;
using std::deque;

class ZergUnit : public BaseUnit
{

public:
//    ZergUnit();

//    ZergUnit(string name);

//    ZergUnit(const UnitStats& baseStats);

//    ZergUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

//    ZergUnit(ZergUnit const& zergUnit);

//    ZergUnit(BaseUnit const& baseUnit);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer <= 0)
        {
            if(mStats.health > EPS && mStats.health < mStats.maxHealth)
            {
                BaseUnit::addHealth(0.27 * state.regenerationUpdate / 1000);
            }
        }
    }
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        regenerate(own);
    }
    //void initUpgrades(vector<int> const& flags);

};

// Zerg Units

class Drone final : public ZergUnit
{};

class Queen final : public ZergUnit
{};

class Zergling final : public ZergUnit
{
public:
    void initUpgrades (vector<int> const& flags);
};

class Baneling final : public ZergUnit
{
public:
    void initUpgrades (vector<int> const& flags);
};

class Roach final : public ZergUnit
{
public:
    void initUpgrades (vector<int> const& flags);
};

class Hydralisk final : public ZergUnit
{};

class Infestor final : public ZergUnit
{};

class SwarmHost final : public ZergUnit
{};

class Ultralisk :public ZergUnit
{};

class Overseer final : public ZergUnit
{};

class Mutalisk final : public ZergUnit
{};

class Corruptor final : public ZergUnit
{};

class BroodLord final : public ZergUnit
{};

class Viper final : public ZergUnit
{};

#endif // ZERG_H
