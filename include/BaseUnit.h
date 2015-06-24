#ifndef BASEUNIT_H_
#define BASEUNIT_H_

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
#include "Chromosome.h"
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

enum Attribute
{
    light = 0, armored = 1, biological = 2, mechanical = 3, psyonic = 4, massive = 5, air = 6
};

struct Bonus
{
    double base;
    double upgrade;
    vector<Attribute> attributes;
    Bonus() : base(0.0), upgrade(0.0)
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
    double minerals = 0;
    double gas = 0;
    double groundRange = 0;
    double airRange = 0;
    double health = 0;
    double maxHealth = 0;
    double shield = 0;
    double maxShield = 0;
    double armor = 0;
    double armorUpgradeBonus = 0;
    double speed = 0;
    double energy = 0;
    double maxEnergy = 0;
    double size = 0;
    bool airUnit = false;

    double groundAttack = 0;
    double gaUpgrade = 0;

    double airAttack = 0;
    double aaUpgrade = 0;

    int gaCooldown = 0;
    int aaCooldown = 0;

    double sumMaxHealthAndShield = 0;

    double creepMultiplier = 1.0;


    vector<Attribute> attributes;
    vector<Bonus> bonuses;

};

class BaseUnit
{  
protected:

    bitset<NBITS> const* mChromosome;
    vector<double> mPhenotype;
    size_t mChromosomeStartPosition;

    UnitStats mStats;

    Vec2D currentForce;
    Vec2D mPos;
    Vec2D mMinPos;
    Vec2D mMaxPos;
    double mMaxDist;
    Vec2D mStartPos;
    bool mTracking = false;
    vector<Vec2Df> mPath;
    vector<BaseUnit *> mInRange;

    int mTimeSlice = 10;

    int mAttackTimer = 0;
    int mMovementTimer = 0;

    int mMovementUpdate = 100;

    int mAttackUpgrade = 0.0;
    int mArmorUpgrade = 0.0;
    int mShieldUpgrade = 0.0;

    double mAttackMultiplier = 1.0;
    double mDefenseMultiplier = 1.0;

    double mDefenseSubtractor = 0.0;


    string mName;

    int mId;

    double mMoveDist;
    double mMovementUpdateDist;

    BaseUnit *mTarget = nullptr;







    std::function<Vec2D(BaseUnit & own, BaseUnit & buddy)> mFriendFunc = [] (BaseUnit & own, BaseUnit & buddy)
    {
        Vec2D distVec(buddy.getX() - own.getX(), buddy.getY() - own.getY());
        double const dist = distVec.computeLength();
        if(dist < own.getSize () + buddy.getSize())
        {
            Vec2D force = std::move(distVec.getNormedVec(dist));
            return Vec2D(-1e5*force.x, -1e5*force.y);
        }


        // tmp[0] = own.getMaxDist()*own.getPhenotype(0)
        if(dist < own.tmp[0])
        {
            Vec2D res = std::move(distVec.getNormedVec(dist));
            // tmp[1] = 1e3*own.getPhenotype(1) + 1e2*own.getResources()*own.getPhenotype(2)
            res.x *= own.tmp[1];
            res.y *= own.tmp[1];
            return res;
        }
        else
        {
            return Vec2D(0.0);
        }
    };



    std::function<Vec2D(BaseUnit & own, BaseUnit & buddy)> mEnemyFunc = [] (BaseUnit & own, BaseUnit & enemy)
    {
        Vec2D distVec(enemy.getX() - own.getX(), enemy.getY() - own.getY());
        double const dist = distVec.computeLength ();
        if(dist < own.getSize() + enemy.getSize())
        {
            Vec2D force = std::move(distVec.getNormedVec(dist));
            return Vec2D(-1e5*force.x, -1e5*force.y);
        }

        Vec2D res1, res2;
        double const ownRange = own.computeRange (enemy);

        int const enemyId = enemy.getIdentifier();
        if(!own.mPossibleDamage[enemyId].valid)
        {
            own.mPossibleDamage[enemyId] = own.computeDamage(enemy);
        }
        Damage const& ownDamage = own.mPossibleDamage[enemyId];
        double const enemyMovement = enemy.getMovementUpdateDist();
        if(dist > ownRange*own.getPhenotype(3))
        {
            res1 = distVec.getNormedVec(dist);
            double const a = enemy.getSumMaxHealthShield () - enemy.getHealth () - enemy.getShield ();
            double const b = enemy.isAirUnit () ? std::max(own.getAACooldown () - own.getAttackTimer(),0)
                                             : std::max(own.getGACooldown () - own.getAttackTimer (),0);
            // tmp[2] = 1e1*own.getPhenotype(4)
            // tmp[3] = 1e2*own.getPhenotype(6)
            // tmp[4] = 1e3*own.getPhenotype(7)
            double const val = a*own.tmp[2] + b*own.getPhenotype(5) + ownDamage.total*own.tmp[3] + own.tmp[4];
            res1.x *= val;
            res1.y *= val;
        }
        else if(dist < (ownRange - enemyMovement)*own.getPhenotype(8))
        {
            res1 = distVec.getNormedVec(dist);
            //tmp[5] = -1e4*own.getPhenotype(9)
            res1.x *= own.tmp[5];
            res1.y *= own.tmp[5];
        }

        double const enemyRange = enemy.computeRange(own);
        int const ownId = own.getIdentifier();
        if(!enemy.mPossibleDamage[ownId].valid)
        {
            enemy.mPossibleDamage[ownId] = enemy.computeDamage(own);
        }
        Damage const& enemyDamage = enemy.mPossibleDamage[ownId];

        if(dist < own.getPhenotype (10)*(enemyRange + enemyMovement))
        {

            res2 = distVec.getNormedVec(dist);
            double const a = own.getSumMaxHealthShield () - own.getHealth () - own.getShield ();
            double const b = own.isAirUnit () ? std::max(enemy.getAACooldown () - enemy.getAttackTimer(),0)
                                           : std::max(enemy.getGACooldown () - enemy.getAttackTimer (),0);
            // tmp[6] = 1e1*own.getPhenotype(11)
            // tmp[7] = 1e2*own.getPhenotype(13)
            // tmp[8] = 1e3*own.getPhenotype(14)
            double const val = a*own.tmp[6] + b*own.getPhenotype(12) + enemyDamage.total*own.tmp[7] + own.tmp[8];
            res2.x *= val;
            res2.y *= val;
        }

        return Vec2D(res1.x-res2.x, res1.y-res2.y);
    };

    Damage computeDamage(BaseUnit const& unit) const;

    bool attack(BaseUnit& unit);

public:

    int const mNGenes = 16;

    Damage mPossibleDamage[18];

    double tmp[9];

    int mMovementUpdateBackup;

    bool mCSAffected = false;

    BaseUnit();

    BaseUnit(string name);

    BaseUnit(const UnitStats& baseStats);

    BaseUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

    //BaseUnit(BaseUnit const& baseUnit);


    Vec2D computeFriendForce(BaseUnit & other);

    Vec2D computeEnemyForce(BaseUnit & other);

    void setName(string name);

    string getName() const;

    // getter and setter
    void setStats(const UnitStats& newStats);

    void setSpeed(double value);

    void setMaxHealth(double const value);

    void incMaxHealth(double const value);

    void incHealth(double const value);

    UnitStats& accessStats();

    double getMaxDist() const;

    double getMinerals() const;

    double getGas() const;

    double getResources() const;

    double getGroundRange() const;

    double getAirRange() const;

    double getHealth() const ;

    double getMaxHealth() const;


    double getShield() const;

    double getMaxShield() const;

    double getSumMaxHealthShield() const;

    double getArmor() const;

    double getSight() const;

    double getSpeed() const;

    double getAcceleration() const;

    double getEnergy() const;

    double getMaxEnergy() const;

    bool isAirUnit() const;

    double getSize() const;

    void addHealth(double const value);

    void subHealth(double const value);

    void setHealth(double const value);

    void addShield(double const value);

    void subShield(double const value);

    void incArmor();

    void decArmor();

    void resetHealth();

    void resetShield();

    void resetEnergy();

    double getGroundAttack() const;

    double getAirAttack() const;

    int getGACooldown() const;

    int getAACooldown() const;

    double getGAUpgradeBonus() const;

    double getAAUpgradeBonus() const;

    double getArmorUpgradeBonus() const;

    int getArmorUpgrade() const;

    int getAttackUpgrade() const;

    vector<Attribute> const& getAttributes() const;

    vector<Bonus> const& getBonuses() const;

    void computeTemporaryValues();



    Vec2D getPos() const;

    Vec2D getMinPos() const;

    Vec2D getMaxPos() const;

    double getX() const;

    double getY() const;

    void setPos(double const x, double const y);

    void setPos(const Vec2D pos);

    void setMinPos(Vec2D const pos);

    void setMaxPos(Vec2D const pos);

    void setPosLimits(Vec2D const minPos, Vec2D const maxPos);

    void setX(double const x);

    void setY(double const y);

    void setFriendForce(function<Vec2D(BaseUnit &, BaseUnit &)> func);

    void setEnemyForce(function<Vec2D(BaseUnit &, BaseUnit &)> func);

    bool isDead() const;

    double computeRange(BaseUnit const& other) const;
    double computeAirRange(BaseUnit const& other) const;
    double computeGroundRange(BaseUnit const& other) const;

    void initUpgrades(vector<int> const& flags);


    template<typename T> double computeDistance(T const& unit)
    {
        double const delta_x = unit.getX()-mPos.x;
        double const delta_y = unit.getY()-mPos.y;
        return std::sqrt(delta_x*delta_x + delta_y*delta_y);
    }

    template<typename T> double computeDistanceSquared(T const& unit)
    {
        double const delta_x = unit.getX()-mPos.x;
        double const delta_y = unit.getY()-mPos.y;
        return delta_x*delta_x + delta_y*delta_y;
    }



    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(!(attack(other)))
        {
            move(own, other);
        }
        decMovementTimer();
        decAttackTimer ();
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
        if(mMovementTimer <= 0)
        {
            Vec2D force(0.0);
            for(const auto& pot : own.potentialList)
            {
                Vec2D const generatedForce = std::move(pot.computeForce(*this));
                force.x += generatedForce.x;
                force.y += generatedForce.y;
            }
            for(auto buddy :  own.unitList)
            {
                if(buddy->isDead() || static_cast<BaseUnit*>(buddy) == this)
                {
                    continue;
                }
                Vec2D const generatedForce = std::move(this->computeFriendForce(*buddy));
                force.x += generatedForce.x;
                force.y += generatedForce.y;
            }
            for(auto enemy :  other.unitList)
            {
                if(enemy->isDead())
                {
                    continue;
                }

                Vec2D const generatedForce = std::move(this->computeEnemyForce(*enemy));
                force.x += generatedForce.x;
                force.y += generatedForce.y;
            }
            for(auto const& forceField : own.forceFieldQueue)
            {
                Vec2D const generatedForce = std::move(forceField.second.computeForce(*this));
                force.x += generatedForce.x;
                force.y += generatedForce.y;
            }

            currentForce = std::move(force.getNormedVec());
            mMovementTimer = mMovementUpdate;
        }

        setPos (currentForce.x*mMoveDist + mPos.x, currentForce.y*mMoveDist + mPos.y);

    }


    bool attackPossible(BaseUnit const& enemy);

    template<typename T> bool attack(PlayerState<T>& other)
    {
        if(this->mAttackTimer > 0 || this->isDead() || other.unitList.empty())
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
            else
            {
                continue;
            }

            bool const newKill = getPhenotype(15) * damage.total > enemy->getHealth() + enemy->getShield();
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

    double getPhenotype(size_t const pos) const;

    void setChromosomeStartPosition(size_t const pos);
    void setChromosome(Chromosome const & chromosome, size_t const pos);
    void setChromosome(Chromosome const & chromosome);

    void setTracking(bool const tracking);
    void reservePathStorage(size_t const sz);
    vector<Vec2Df> getPath() const;
    void clearPath();

    int getAttackTimer() const;
    void setAttackTimer(int value);
    void decAttackTimer();


    void resetTimer();

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

    Vec2D getStartPos() const;
    void setStartPos(const Vec2D &value);
    void resetPos();
    int getMovementTimer() const;
    void setMovementTimer(int value);
    void decMovementTimer();
    int getMovementUpdate() const;
    void setMovementUpdate(int value);

    double getMoveDist() const;
    double getMovementUpdateDist() const;
    void multSpeed(double value);

    size_t getChromosomeStartPosition() const;
};


/* Race Base Unit Types, as there exist unique unit characteristics for the three different races
 f.e. Zerg Units regenerate health and Protoss Units have a regenerating shield */



#endif


