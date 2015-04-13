#ifndef _UNIT_H_
#define _UNIT_H_

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


#include "PlayerState.h"
#include "Utilities.h"
#include "UnitGenes.h"

using std::pair;
using std::tuple;
using std::function;
using std::pow;
using std::sqrt;
using std::list;
using std::string;
using std::vector;
using std::array;


enum Attribute
{
    light = 0, armored = 1, biological = 2, mechanical = 3, psyonic = 4, massive = 5, air = 6
};

struct Bonus
{
    float base;
    float upgrade;
    vector<Attribute> attributes;
    Bonus() : base(0.f), upgrade(0.f)
    {}
    Bonus(Bonus const& other) : base(other.base), upgrade(other.upgrade), attributes(other.attributes)
    {}
};

struct Damage
{
    double shield;
    double health;
    double total;
    bool valid;

    Damage() : shield(0.0), health(0.0), total(0.0), valid(false) {}
    Damage(double const s, double const h, double const t) : shield(s), health(h), total(t), valid(true) {}
    Damage(Damage const& damage) : shield(damage.shield), health(damage.health), total(damage.total), valid(damage.valid) {}
};


struct UnitStats
{
    float minerals = 0;
    float gas = 0;
    float gdps = 0;
    float adps = 0;
    float groundRange = 0;
    float airRange = 0;
    double health = 0;
    double maxHealth = 0;
    double shield = 0;
    double maxShield = 0;
    float armor = 0;
    float armorUpgradeBonus = 0;
    float speed = 0;
    double energy = 0;
    double maxEnergy = 0;
    float size = 0;
    bool airUnit = false;

    float groundAttack = 0;
    float gaUpgradeBonus = 0;

    float airAttack = 0;
    float aaUpgradeBonus = 0;
    int gaCooldown = 0;
    int aaCooldown = 0;


    vector<Attribute> attributes;
    vector<Bonus> bonuses;

};

class BaseUnit
{  
protected:

    UnitGenes mGenes;
    UnitStats mStats;

    Vec2D currentForce;
    Vec2D mPos;
    Vec2D mMinPos;
    Vec2D mMaxPos;
    bool mTracking = false;
    vector<Vec2Df> mPath;
    vector<BaseUnit *> mInRange;

    int mTimeSlice = 10;

    int mAttackTimer = 0;

    int mAttackUpgrade = 0;
    int mArmorUpgrade = 0;
    int mShieldUpgrade = 0;

    double mAttackMultiplier = 1.0;
    double mDefenseMultiplier = 1.0;

    double mDefenseSubtractor = 0.0;

    string mName;

    int mId;







    std::function<Vec2D(BaseUnit & own, BaseUnit & buddy)> mFriendFunc = [] (BaseUnit & own, BaseUnit & buddy)
    {
       Vec2D distVec(buddy.getX() - own.getX(), buddy.getY() - own.getY());

       if(distVec.computeLength() < own.getMaxDist()*(static_cast<double>(own.getYGene(0))/YMAX))
       {
           Vec2D res = distVec.getNormedVec();
           double const val = own.getXGene(0)*buddy.getResources()+own.getYGene(1);
           res.x *= val;
           res.y *= val;
           return res;
       }
       else
       {
           return Vec2D(0.);
       }
    };



    std::function<Vec2D(BaseUnit & own, BaseUnit & buddy)> mEnemyFunc = [] (BaseUnit & own, BaseUnit & enemy)
    {

        Vec2D distVec(enemy.getX() - own.getX(), enemy.getY() - own.getY());
        double const ownRange = own.computeRange(enemy);
        Vec2D res1;

        int const ownId = own.getIdentifier();
        int const enemyId = enemy.getIdentifier();
        if(!own.mPossibleDamage[enemyId].valid)
        {
            own.mPossibleDamage[enemyId] = own.computeDamage(enemy);
        }
        if(!enemy.mPossibleDamage[ownId].valid)
        {
            enemy.mPossibleDamage[ownId] = enemy.computeDamage(own);
        }
        Damage const& enemyDamage = enemy.mPossibleDamage[ownId];
        Damage const& ownDamage = own.mPossibleDamage[enemyId];

        double const len = distVec.computeLength ();
        if(len > ownRange)
        {
            res1 = distVec.getNormedVec();
            double const val = own.getXGene(1)*enemy.getResources()+own.getXGene(2)*ownDamage.total+own.getYGene(2);
            res1.x *= val;
            res1.y *= val;
        }
        else if(len < ownRange-(ownRange*static_cast<double>(own.getYGene(3))/YMAX/4.))
        {
            res1 = distVec.getNormedVec();
            res1.x *= -own.getYGene(4);
            res1.y *= -own.getYGene(4);
        }
        else
        {
            res1 = Vec2D(0.0);
        }
        Vec2D res2;
        double const enemyRange = enemy.computeRange(own);
        if(len < enemyRange+enemyRange*static_cast<double>(own.getYGene(5))/YMAX/4.)
        {
            res2 = distVec.getNormedVec();
            double const val = -own.getYGene(6)-own.getXGene(3)*own.getResources()-own.getXGene(4)*enemyDamage.total;
            res2.x *= val;
            res2.y *= val;
        }
        else
        {
            res2 = Vec2D(0.0);
        }

        return Vec2D(res1.x+res2.x, res1.y+res2.y);
    };

    Damage computeDamage(BaseUnit const& unit) const;

    bool attack(BaseUnit& unit);

public:

    Damage mPossibleDamage[18];

    BaseUnit();

    BaseUnit(string name);

    BaseUnit(const UnitStats& baseStats);

    BaseUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

    BaseUnit(BaseUnit const& baseUnit);


    Vec2D computeFriendForce(BaseUnit & other);

    Vec2D computeEnemyForce(BaseUnit & other);

    void setName(string name);

    string getName() const;

    // getter and setter
    void setStats(const UnitStats& newStats);

    double getMaxDist() const;

    float getMinerals() const;

    float getGas() const;

    float getResources() const;

    float getGdps() const ;

    float getAdps() const;

    float getGroundRange() const;

    float getAirRange() const;

    double getHealth() const ;

    double getMaxHealth() const;

    double getShield() const;

    double getMaxShield() const;

    float getArmor() const;

    float getSight() const;

    float getSpeed() const;

    float getAcceleration() const;

    double getEnergy() const;

    double getMaxEnergy() const;

    bool isAirUnit() const;

    float getSize() const;

    void addHealth(double const value);

    void subHealth(double const value);

    void addShield(double const value);

    void subShield(double const value);

    void incArmor();

    void decArmor();

    void resetHealth();

    void resetShield();

    void resetEnergy();

    float getGroundAttack() const;

    float getAirAttack() const;

    int getGACooldown() const;

    int getAACooldown() const;

    float getGAUpgradeBonus() const;

    float getAAUpgradeBonus() const;

    float getArmorUpgradeBonus() const;

    int getArmorUpgrade() const;

    int getAttackUpgrade() const;

    vector<Attribute> const& getAttributes() const;

    vector<Bonus> const& getBonuses() const;



    Vec2D getPos() const;

    Vec2D getMinPos() const;

    Vec2D getMaxPos() const;

    double getX() const;

    double getY() const;

    void setPos(double const x, double const y);

    void setPos(const Vec2D pos);

    void setMinPos(Vec2D const pos);

    void setMaxPos(Vec2D const pos);

    void setX(double const x);

    void setY(double const y);

    void setFriendForce(function<Vec2D(BaseUnit &, BaseUnit &)> func);

    void setEnemyForce(function<Vec2D(BaseUnit &, BaseUnit &)> func);

    double computeRange(BaseUnit const& other) const;



    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(!(attack(other)))
        {
            move(own, other);
        }
        decAttackTimer();
    }



    template <typename T, typename U> void move(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(this->getHealth () < EPS)
        {
            return;
        }
        if(mTracking)
        {
            mPath.emplace_back(static_cast<float>(mPos.x),static_cast<float>(mPos.y));
        }
        if(own.movementTimer <= 0)
        {
            Vec2D force(0.0);
            for(const auto& pot : own.potentialList)
            {
                const Vec2D generatedForce = pot.computeForce(*this);
                force.x += generatedForce.x;
                force.y += generatedForce.y;
            }
            for(auto buddy :  own.unitList)
            {
                if(buddy->getHealth() < EPS)
                {
                    continue;
                }

                const Vec2D generatedForce = this->computeFriendForce(*buddy);
                force.x += generatedForce.x;
                force.y += generatedForce.y;
            }
            for(auto enemy :  other.unitList)
            {
                if(enemy->getHealth() < EPS)
                {
                    continue;
                }

                const Vec2D generatedForce = this->computeEnemyForce(*enemy);
                force.x += generatedForce.x;
                force.y += generatedForce.y;
            }
            currentForce = force.getNormedVec ();
        }
        double const moveDist = getSpeed()*own.movementUpdate/1000;


        setX(currentForce.x*moveDist+getX());
        setY(currentForce.y*moveDist+getY());

        // Collision detection, include it again if needed
        //TODO needs to be fixed before included again!
        /*
        for(auto unit : other.unitList)
        {
            double x = this->getX()-unit->getX();
            double y = this->getY()-unit->getY();
            double const dist = std::sqrt(pow(x,2)+pow(y,2));
            double const allowed_dist = unit->getSize() + this->getSize();
            if(dist < allowed_dist)
            {
                double len = allowed_dist - dist;
                Vec2D shift = normVecSafe(Vec2D(x, y), dist);
                setX(getX()+len*shift.x);
                setY(getY()+len*shift.y);
            }

        }
        for(auto unit : own.unitList)
        {
            double x = this->getX()-unit->getX();
            double y = this->getY()-unit->getY();
            double const dist = std::sqrt(pow(x,2)+pow(y,2));
            double const allowed_dist = unit->getSize() + this->getSize();
            if(dist < allowed_dist)
            {
                double len = allowed_dist - dist;
                Vec2D shift = normVecSafe(Vec2D(x, y), dist);
                setX(getX()+len*shift.x);
                setY(getY()+len*shift.y);
            }

        }
        */

    }

    bool attackPossible(BaseUnit const& enemy);

    template<typename T> bool attack(PlayerState<T>& other)
    {
        if(this->mAttackTimer > 0 || this->getHealth() < EPS || other.unitList.empty())
        {
            return false;
        }

        auto aim = *(other.unitList.begin());
        double maxDamage = 0.0;
        bool kill = false;
        for(auto enemy : other.unitList)
        {
            Damage damage;
            if(attackPossible (*enemy))
            {
                damage = computeDamage(*enemy);
            }

            bool const newKill = 0.75 * damage.total > enemy->getHealth() + enemy->getShield();
            bool const higher = damage.total > maxDamage;
            if(newKill)
            {
                if(kill)
                {
                    if(higher)
                    {
                        maxDamage = damage.total;
                        aim = enemy;
                        mPossibleDamage[enemy->getIdentifier()] = damage;
                    }
                }
                else
                {
                    maxDamage = damage.total;
                    aim = enemy;
                    mPossibleDamage[enemy->getIdentifier()] = damage;
                }
            }
            else
            {
                if(higher)
                {
                    maxDamage = damage.total;
                    aim = enemy;
                    mPossibleDamage[enemy->getIdentifier()] = damage;
                }
            }
        }
        if(maxDamage < EPS)
        {
            return false;
        }
        if(attack(*aim))
        {
            --other.unitCount;
        }
        return true;
    }

    int getXGene(size_t pos) const;
    int getYGene(size_t pos) const;
    void setGenes(UnitGenes const& genes);
    size_t getHash() const;

    void setTracking(bool const tracking);
    void reservePathStorage(size_t const sz);
    vector<Vec2Df> getPath() const;
    void clearPath();

    int getAttackTimer() const;
    void setAttackTimer(int value);
    void decAttackTimer();

    int getTimeSlice() const;
    void setTimeSlice(int value);

    int getShieldUpgrade() const;
    void setShieldUpgrade(int value);
    double getAttackMultiplier() const;
    void setAttackMultiplier(double value);
    double getDefenseMultiplier() const;
    void setDefenseMultiplier(double value);
    double getDefenseSubtractor() const;
    void setDefenseSubtractor(double value);
    int getIdentifier() const;
    void setIdentifier(int value);
};


/* Race Base Unit Types, as there exist unique unit characteristics for the three different races
 f.e. Zerg Units regenerate health and Protoss Units have a regenerating shield */
class TerranUnit : public BaseUnit
{
public:
    TerranUnit();

    TerranUnit(string name);

    TerranUnit(const UnitStats& baseStats);

    TerranUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

    TerranUnit(TerranUnit const& terranUnit);

    TerranUnit(BaseUnit const& baseUnit);

    template <typename T> void regenerate(PlayerState<T>&) {}

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        regenerate(own);
    }
};

class ZergUnit : public BaseUnit
{

public:
    ZergUnit();

    ZergUnit(string name);

    ZergUnit(const UnitStats& baseStats);

    ZergUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

    ZergUnit(ZergUnit const& zergUnit);

    ZergUnit(BaseUnit const& baseUnit);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer == 0)
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

};

class ProtossUnit : public BaseUnit
{
protected:
    int mShieldRegenCount = 0;
public:
    ProtossUnit();

    ProtossUnit(string name);

    ProtossUnit(const UnitStats& baseStats);

    ProtossUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

    ProtossUnit(ProtossUnit const& protossUnit);

    ProtossUnit(BaseUnit const& baseUnit);

    void subShield(double const value);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer == 0)
        {
            if(mShieldRegenCount == 0 && mStats.shield < mStats.maxShield && mStats.health > EPS)
            {
                BaseUnit::addShield(2.0 * state.regenerationUpdate / 1000);
            }
            else
            {
                --mShieldRegenCount;
            }
        }

    }

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        BaseUnit::timestep(own, other);
        regenerate(own);
    }
};



// classes for every specific Unit

// Terran Units

class SCV final : public TerranUnit
{};

class Marine final : public TerranUnit
{};

class Marauder final : public TerranUnit
{};

class Reaper final : public TerranUnit
{
protected:
    int mHealthRegenCount = 0;
public:

    Reaper();

    Reaper(string name);

    Reaper(const UnitStats& baseStats);

    Reaper(const UnitStats& baseStats, Vec2D min, Vec2D max);

    Reaper(Reaper const& Reaper);

    Reaper(BaseUnit const& baseUnit);

    Reaper(TerranUnit const& terranUnit);

    void subHealth(double const value);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer == 0)
        {
            if(mHealthRegenCount == 0 && mStats.health < mStats.maxHealth && mStats.health > EPS)
            {
                BaseUnit::addHealth(2.0 * state.regenerationUpdate / 1000);
            }
            else
            {
                --mHealthRegenCount;
            }
        }
    }

    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        TerranUnit::timestep(own, other);
        regenerate(own);
    }
};

class Ghost final : public TerranUnit
{};

class Hellion final : public TerranUnit
{};

class Hellbat final : public TerranUnit
{};

class SiegeTank final : public TerranUnit
{};

class WidowMine final : public TerranUnit
{};

class Thor final : public TerranUnit
{};

class Viking final : public TerranUnit
{};

class Medivac final : public TerranUnit
{};

class Raven final : public TerranUnit
{};

class Banshee final : public TerranUnit
{};

class Battlecruiser final : public TerranUnit
{};

// Zerg Units

class Drone final : public ZergUnit
{};

class Queen final : public ZergUnit
{};

class Zergling final : public ZergUnit
{};

class Baneling final : public ZergUnit
{};

class Roach final : public ZergUnit
{};

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


// Protoss Units

class Probe final : public ProtossUnit
{};

class Zealot final : public ProtossUnit
{};

class Stalker final : public ProtossUnit
{};

class Sentry final : public ProtossUnit
{};

class HighTemplar final : public ProtossUnit
{};

class DarkTemplar final : public ProtossUnit
{};

class Immortal final : public ProtossUnit
{};

class Colossus final : public ProtossUnit
{};

class Archon final : public ProtossUnit
{};

class Observer final : public ProtossUnit
{};

class WarpPrism final : public ProtossUnit
{};

class Phoenix final : public ProtossUnit
{};

class VoidRay final : public ProtossUnit
{};

class Oracle final : public ProtossUnit
{};

class Carrier final : public ProtossUnit
{};

class Tempest final : public ProtossUnit
{};

class MothershipCore final : public ProtossUnit
{};

class Mothership final : public ProtossUnit
{};


#endif


