#ifndef BASEUNIT_H_
#define BASEUNIT_H_

#include<utility>
#include<functional>
#include<cmath>
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


// Unit attributes
enum Attribute
{
    light = 0, armored = 1, biological = 2, mechanical = 3, psyonic = 4, massive = 5, air = 6
};

struct Bonus
{
    double base;
    double upgrade;
    std::vector<Attribute> attributes;
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

    std::vector<Attribute> attributes;
    std::vector<Bonus> bonuses;
};


// Basic Unit class that contains all shared functionality

class BaseUnit
{  
protected:

    // Pointer to the combined chromosome of all units (stored in PlayerState)
    double const* mChromosome;

    // Start position of the units actual chromosome in the combined one
    size_t mChromosomeStartPosition;

    UnitStats mStats;

    Vec2D currentForce;
    Vec2D mPos;
    Vec2D mMinPos;
    Vec2D mMaxPos;
    double mMaxDist;
    Vec2D mStartPos;

    // Should the path of the unit be saved
    bool mTracking = false;
    std::vector<Vec2Df> mPath;

    std::vector<BaseUnit *> mInRange;

    // Length of the time step
    int mTimeSlice = 10;

    // Timer for ensuring that the unit only repeatedly attacks after the cooldown hast passed
    int mAttackTimer = 0;

    // Timer for recomputing the movement direction
    int mMovementTimer = 0;

    int mMovementUpdate = 100;

    int mAttackUpgrade = 0.0;
    int mArmorUpgrade = 0.0;
    int mShieldUpgrade = 0.0;

    double mAttackMultiplier = 1.0;
    double mDefenseMultiplier = 1.0;

    double mDefenseSubtractor = 0.0;


    std::string mName;

    // Identifier for fast identification of specific units
    int mId;

    double mMoveDist;
    double mMovementUpdateDist;

    // Last attack target of the unit
    BaseUnit *mTarget = nullptr;

    bool mHasAttacked = false;

    bool mMelee = false;


//    tmp[0] = getMaxDist()*getPhenotype(1);
//    tmp[1] = getMaxDist()*getPhenotype(2) + getSize();
//    tmp[2] = 1e3*getPhenotype(3) + 1e1*getResources()*getPhenotype(4);
//    tmp[3] = 1e1*getPhenotype(6);
//    tmp[4] = 1e2*getPhenotype(8);
//    tmp[5] = 1e3*getPhenotype(9);
//    tmp[6] = 1e1*getPhenotype(12);
//    tmp[7] = 1e2*getPhenotype(13);
//    tmp[8] = 1e3*getPhenotype(14);


    // Computes the force resulting from the attractive potential of friendly units

    std::function<Vec2D(BaseUnit & own, BaseUnit & buddy)> mFriendFunc = [] (BaseUnit & own, BaseUnit & buddy)
    {
        Vec2D distVec(buddy.getX() - own.getX(), buddy.getY() - own.getY());
        double const dist = distVec.computeLength();

        // Collision avoidance
        if(dist < own.getSize () + buddy.getSize())
        {
            Vec2D force = std::move(distVec.getNormedVec(dist));
            return Vec2D(-1e6*force.x, -1e6*force.y);
        }

        if(dist < own.tmp[0] && dist > own.tmp[1] + buddy.getSize())
        {
            Vec2D res = std::move(distVec.getNormedVec(dist));
            res.x *= own.tmp[2];
            res.y *= own.tmp[2];
            return res;
        }
        else
        {
            return Vec2D(0.0);
        }
    };


    // Computes the force resulting from the attractive and repulsive potential of enemy units

    std::function<Vec2D(BaseUnit & own, BaseUnit & buddy)> mEnemyFunc = [] (BaseUnit & own, BaseUnit & enemy)
    {
        Vec2D distVec(enemy.getX() - own.getX(), enemy.getY() - own.getY());
        double const dist = distVec.computeLength ();
        // Collision avoidance
        if(dist < own.getSize() + enemy.getSize())
        {
            Vec2D force = std::move(distVec.getNormedVec(dist));
            return Vec2D(-1e6*force.x, -1e6*force.y);
        }

        Vec2D normedDistVec = std::move(distVec.getNormedVec(dist));

        Vec2D res1;
        double const ownRange = own.computeRange (enemy);

        // Attraction
        if(dist > ownRange*own.getPhenotype(5))
        {
            int const enemyId = enemy.getIdentifier();
            if(!own.mPossibleDamage[enemyId].valid)
            {
                own.mPossibleDamage[enemyId] = own.computeDamage(enemy);
            }
            Damage const& ownDamage = own.mPossibleDamage[enemyId];
            double const a = enemy.getSumMaxHealthShield () - enemy.getHealth () - enemy.getShield ();
            double const b = enemy.isAirUnit () ? std::max(own.getAACooldown () - own.getAttackTimer(),0)
                                             : std::max(own.getGACooldown () - own.getAttackTimer (),0);
            res1 = normedDistVec;
            double const val = a*own.tmp[3] + b * own.getPhenotype(7) + ownDamage.total*own.tmp[4] + own.tmp[5];
            res1.x *= val;
            res1.y *= val;
        }

        Vec2D res2;
        double const enemyRange = enemy.computeRange(own);
        double const enemyMovement = enemy.getMovementUpdateDist();

        // Repulsion
        if(dist < own.getPhenotype (10)*(enemyRange+enemyMovement))
        {
            res2 = normedDistVec;
            int const ownId = own.getIdentifier();
            if(!enemy.mPossibleDamage[ownId].valid)
            {
                enemy.mPossibleDamage[ownId] = enemy.computeDamage(own);
            }
            Damage const& enemyDamage = enemy.mPossibleDamage[ownId];

            double const a = own.getSumMaxHealthShield () - own.getHealth () - own.getShield ();
            double const b = own.isAirUnit () ? std::max(enemy.getAACooldown () - enemy.getAttackTimer(),0)
                                           : std::max(enemy.getGACooldown () - enemy.getAttackTimer (),0);


            double const val = a*own.tmp[6] + b*own.getPhenotype(11) + enemyDamage.total*own.tmp[7] + own.tmp[8];
            res2.x *= val;
            res2.y *= val;
        }


        return Vec2D(res1.x-res2.x, res1.y-res2.y);
    };

    // Compute the damage the unit could theoretically cause to a unit
    Damage computeDamage(BaseUnit const& unit) const;

    // Attack a specific unit
    bool attack(BaseUnit& unit);

public:

    // Length of the chromosome
    int const mNGenes = 15;

    // Possible damage that can be applied to units with identifier = position in array
    Damage mPossibleDamage[18];

    // Intermediate results of the potenial field computation
    double tmp[10];

    int mMovementUpdateBackup;

    // Is the unit affected by Concussive Shells (Marauder ability)
    bool mCSAffected = false;

    BaseUnit();

    BaseUnit(std::string name);

    BaseUnit(const UnitStats& baseStats);

    BaseUnit(const UnitStats& baseStats, Vec2D min, Vec2D max);

    //BaseUnit(BaseUnit const& baseUnit);


    Vec2D computeFriendForce(BaseUnit & other);

    Vec2D computeEnemyForce(BaseUnit & other);

    void setName(std::string name);

    std::string getName() const;

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

    std::vector<Attribute> const& getAttributes() const;

    std::vector<Bonus> const& getBonuses() const;

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

    void setFriendForce(std::function<Vec2D(BaseUnit &, BaseUnit &)> func);

    void setEnemyForce(std::function<Vec2D(BaseUnit &, BaseUnit &)> func);

    bool isDead() const;

    double computeRange(BaseUnit const& other) const;
    double computeAirRange(BaseUnit const& other) const;
    double computeGroundRange(BaseUnit const& other) const;

    void initUpgrades(std::vector<int> const& flags);


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



    // Main method used for executing the actions of the next time step
    // Overloaded in subclasses that represent specific units
    template <typename T, typename U> void timestep(PlayerState<T>& own, PlayerState<U>& other)
    {
        if(!(attack(other)))
        {
            move(own, other);
        }
        decMovementTimer();
        decAttackTimer ();
    }


    // Method that emulates the movement of a unit
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
            if(mMelee)
            {
                Vec2D force(0.0);
                double max_diff = std::numeric_limits<double>::min();
                double min_dist = std::numeric_limits<double>::max();
                double min_speed = std::numeric_limits<double>::max();
                bool collision = false;

                for(const auto& pot : own.potentialList)
                {
                    Vec2D const generatedForce = std::move(pot.computeForce(*this));
                    if(std::abs(generatedForce.x) > 1e5 || std::abs(generatedForce.y) > 1e5)
                    {
                        collision = true;
                        force.x += generatedForce.x;
                        force.y += generatedForce.y;
                    }
                }
                for(auto buddyPtr :  own.unitList)
                {
                    auto& buddy = *buddyPtr;
                    if(buddy.isDead() || static_cast<BaseUnit*>(buddyPtr) == this)
                    {
                        continue;
                    }
                    Vec2D distVec(buddy.getX() - getX(), buddy.getY() - getY());
                    double const dist = distVec.computeLength ();
                    // Collision avoidance
                    if(dist < getSize() + buddy.getSize())
                    {
                        collision = true;
                        Vec2D tmp = std::move(distVec.getNormedVec(dist));
                        tmp = std::move(Vec2D(-1e6*tmp.x, -1e6*tmp.y));
                        force.x += tmp.x;
                        force.y += tmp.y;
                    }
                }
                for(auto enemyPtr : other.unitList)
                {
                    auto& enemy = *enemyPtr;
                    if(enemy.isDead())
                    {
                        continue;
                    }
                    Vec2D distVec(enemy.getX() - getX(), enemy.getY() - getY());
                    double const dist = distVec.computeLength ();
                    // Collision avoidance
                    if(dist < getSize() + enemy.getSize())
                    {
                        collision = true;
                        Vec2D tmp = std::move(distVec.getNormedVec(dist));
                        tmp = std::move(Vec2D(-1e6*tmp.x, -1e6*tmp.y));
                        force.x += tmp.x;
                        force.y += tmp.y;
                    }
                    else if(!collision && (enemy.getSpeed() < this->getSpeed() || enemy.getSpeed() < min_speed))
                    {
                        Damage dmg_own = this->computeDamage(enemy);
                        Damage dmg_enemy = enemy.computeDamage(*this);
                        double const diff = dmg_own.total - dmg_enemy.total;
                        if(dist < min_dist || (std::abs(dist - min_dist) < EPS && max_diff < diff))
                        {
                            force = distVec;
                            max_diff = diff;
                            min_speed = enemy.getSpeed();
                            min_dist = dist;
                        }
                    }
                }
                currentForce = std::move(force.getNormedVec());
                mMovementTimer = mMovementUpdate;
            }
            else
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
        }

        setPos (currentForce.x*mMoveDist + mPos.x, currentForce.y*mMoveDist + mPos.y);

    }


    bool attackPossible(BaseUnit const& enemy);

    // Method that emulates that the unit attacks
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

            bool const newKill = getPhenotype(0) * damage.total > enemy->getHealth() + enemy->getShield();
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

    // Methods for manipulating the chromosome of the unit
    double getPhenotype(size_t const pos) const;
    size_t getChromosomeStartPosition() const;
    void setChromosomeStartPosition(size_t const pos);
    void setChromosome(Chromosome const & chromosome, size_t const pos);
    void setChromosome(Chromosome const & chromosome);

    // Methods for enabling the tracking of unit paths
    // currently not used
    void setTracking(bool const tracking);
    void reservePathStorage(size_t const sz);
    std::vector<Vec2Df> getPath() const;
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

    // Change the speed of a unit by a certain factor
    void multSpeed(double value);

    // Has the unit attacked at least once
    bool hasAttacked() const;

    // Reset the unit to the initial state
    void reset();
};


/* Race Base Unit Types, as there exist unique unit characteristics for the three different races
 f.e. Zerg Units regenerate health and Protoss Units have a regenerating shield */



#endif


