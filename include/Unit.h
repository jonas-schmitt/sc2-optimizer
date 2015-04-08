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
    Damage() : shield(0.0), health(0.0), total(0.0) {}
    Damage(double const s, double const h, double const t) : shield(s), health(h), total(t) {}
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

    pair<double,double> currentForce;
    pair<double,double> mPos;
    pair<double,double> mMinPos;
    pair<double,double> mMaxPos;
    bool mTracking = false;
    vector<pair<float,float>> mPath;
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

    std::function<double(BaseUnit const& own, BaseUnit const& buddy, std::vector<std::vector<double>> &field, double off)> mFriendFuncPot = [] (BaseUnit const& own, BaseUnit const& buddy, std::vector<std::vector<double>> &field, double off)
    {
       double const x = buddy.getX()-own.getX();
       double const y = buddy.getY()-own.getY();
       double  dist = std::sqrt(pow(x,2)+pow(y,2));
       pair<double,double> const minPos = buddy.getMinPos();
       pair<double,double> const maxPos = buddy.getMaxPos();

	   double val;
       if(dist < own.getMaxDist()*(static_cast<double>(own.getYGene(0))/YMAX))
       {
           val = own.getXGene(0)*buddy.getResources()+own.getYGene(1);
	   		//val *= (dist+own.getMaxDist()*own.getYGene(0)/YMAX);
		   val *= dist;
	   	//val *= (dist-own.getMaxDist()*own.getYGene(0)/YMAX);
           for (int i = static_cast<int>(std::round(std::max(minPos.first, buddy.getPos().first - buddy.getSize()))); i < static_cast<int>(std::round(std::min(maxPos.first, buddy.getPos().first+buddy.getSize()))); ++i)
		   {
                for (int j = static_cast<int>(std::round(std::max(minPos.second, buddy.getPos().second - buddy.getSize()))); j < static_cast<int>(std::round(std::min(maxPos.second, buddy.getPos().second+buddy.getSize()))); ++j)
				{
					field.at(i).at(j) += (val/off);
				}
				//v.at(buddy->getPos().first).at(buddy->getPos().second) += mFriendFuncPot(*this,*buddy);
			}
	   }
       else
       {
		   val = 0.0f;
       }
	   //val *= (dist-own.getYGene(0));
	   return val;
    };

    std::function<pair<double,double>(BaseUnit const& own, BaseUnit const& buddy)> mFriendFunc = [] (BaseUnit const& own, BaseUnit const& buddy)
    {
       double const x = buddy.getX()-own.getX();
       double const y = buddy.getY()-own.getY();
       double const dist =  std::sqrt(pow(x,2)+pow(y,2));

       /*if (dist < (own.getSize()+buddy.getSize())*1.0f)
	   {
           pair<double,double> res = own.normVecSafe(pair<double,double>(x,y), dist);
		   double const pot = -LIMIT;
		   res.first *= pot;
		   res.second *= pot;
		   return res;
       }*/

       if(dist < own.getMaxDist()*(static_cast<double>(own.getYGene(0))/YMAX))
       {
           pair<double,double> res = own.normVecSafe(pair<double,double>(x,y), dist);
           double const val = own.getXGene(0)*buddy.getResources()+own.getYGene(1);
           res.first *= val;
           res.second *= val;
           return res;
       }
       else
       {
           return pair<double,double>(0.,0.);
       }
    };

    std::function<double(BaseUnit const& own, BaseUnit const& enemy, std::vector<std::vector<double>> &field, double off)> mEnemyFuncPot = [] (BaseUnit const& own, BaseUnit const& enemy, std::vector<std::vector<double>> &field, double off)
	{
        double const x = enemy.getX()-own.getX();
        double const y = enemy.getY()-own.getY();
        double dist =  std::sqrt(pow(x,2)+pow(y,2));
        double const ownRange = own.computeRange(enemy);
       	double val=0.0f;
		double val2 = 0.0f;
        pair<double,double> const minPos = enemy.getMinPos();
        pair<double,double> const maxPos = enemy.getMaxPos();
        Damage enemyDamage = enemy.computeTheoreticalDamage(own);
        Damage ownDamage = own.computeTheoreticalDamage(enemy);
        if(dist > ownRange)
        {
            val = own.getXGene(1)*enemy.getResources()+own.getXGene(2)*ownDamage.total+own.getYGene(2);
			val *= (dist-ownRange);
        }
        else if(dist < ownRange-(ownRange*static_cast<double>(own.getYGene(3))/YMAX/4.))
        {
            val = -own.getYGene(4);
			//val *= (dist-own.getYGene(3)/YMAX/4.);
			val *= dist;
        }
        else
        {
			val = 0.0f;
        }
		//val *= dist;
           for (int i = static_cast<int>(std::round(std::max(minPos.first, enemy.getPos().first - enemy.getSize()))); i < static_cast<int>(std::round(std::min(maxPos.first, enemy.getPos().first+enemy.getSize()))); ++i)
		   {
                for (int j = static_cast<int>(std::round(std::max(minPos.second, enemy.getPos().second - enemy.getSize()))); j < static_cast<int>(std::round(std::min(maxPos.second, enemy.getPos().second+enemy.getSize()))); ++j)
                {
					field.at(i).at(j) += (val/off);
				}
				//v.at(buddy->getPos().first).at(buddy->getPos().second) += mFriendFuncPot(*this,*buddy);
			}
        double const enemyRange = enemy.computeRange(own);
        if(dist < enemyRange+enemyRange*static_cast<double>(own.getYGene(5))/YMAX/4.)
        {
            val2 = (-own.getYGene(6)-own.getXGene(3)*own.getResources()-own.getXGene(4)*enemyDamage.total);
			//val2 *= dist;
			val2 *= (dist-(enemyRange+own.getYGene(5)/YMAX/4.));
           for (int i = static_cast<int>(std::round(std::max(minPos.first,  enemy.getPos().first - enemy.getGroundRange()))); i < static_cast<int>(std::round(std::min(maxPos.first, enemy.getPos().first+enemy.getGroundRange()))); ++i)
		   {
                for (int j = static_cast<int>(std::round(std::max(minPos.second,  enemy.getPos().second - enemy.getGroundRange()))); j < static_cast<int>(std::round(std::min(maxPos.second, enemy.getPos().second+enemy.getGroundRange()))); ++j)
				{
					field.at(i).at(j) += (val2/off);
				}
				//v.at(buddy->getPos().first).at(buddy->getPos().second) += mFriendFuncPot(*this,*buddy);
			}
        }
        else
        {
        }
		//val2 /= 1e4;
		//val /= 1e4;
		return std::max(-1e10, std::min(1e10,val+val2));

	};

    std::function<pair<double,double>(BaseUnit const& own, BaseUnit const& buddy)> mEnemyFunc = [] (BaseUnit const& own, BaseUnit const& enemy)
    {
        double const x = enemy.getX()-own.getX();
        double const y = enemy.getY()-own.getY();
        double const dist =  std::sqrt(pow(x,2)+pow(y,2));
        double const ownRange = own.computeRange(enemy);
        pair<double,double> res1;
        /*
	   if (dist < (own.getSize()+enemy.getSize())*1.0f)
	   {
           pair<double,double> res = own.normVecSafe(pair<double,double>(x,y), dist);
		   double const pot = -LIMIT;
		   res.first *= pot;
		   res.second *= pot;
		   return res;
       }*/
        Damage enemyDamage = enemy.computeTheoreticalDamage(own);
        Damage ownDamage = own.computeTheoreticalDamage(enemy);
        if(dist > ownRange)
        {
            res1 = own.normVecSafe(pair<double,double>(x,y), dist);
            double const val = own.getXGene(1)*enemy.getResources()+own.getXGene(2)*ownDamage.total+own.getYGene(2);
            res1.first *= val;
            res1.second *= val;
        }
        else if(dist < ownRange-(ownRange*static_cast<double>(own.getYGene(3))/YMAX/4.))
        {
            res1 = own.normVecSafe(pair<double,double>(x,y), dist);
            res1.first *= -own.getYGene(4);
            res1.second *= -own.getYGene(4);
        }
        else
        {
            res1 = pair<double,double>(0.,0.);
        }
        pair<double,double> res2;
        double const enemyRange = enemy.computeRange(own);
        if(dist < enemyRange+enemyRange*static_cast<double>(own.getYGene(5))/YMAX/4.)
        {
            res2 = own.normVecSafe(pair<double,double>(x,y), dist);
            double const val = -own.getYGene(6)-own.getXGene(3)*own.getResources()-own.getXGene(4)*enemyDamage.total;
            res2.first *= val;
            res2.second *= val;
        }
        else
        {
            res2 = pair<double,double>(0.,0.);
        }

        return pair<double,double>(res1.first+res2.first, res1.second+res2.second);
    };

    Damage computeTheoreticalDamage(BaseUnit const& unit) const;
    Damage computeDamage(BaseUnit const& unit) const;

    bool attack(BaseUnit& unit);

public:



    BaseUnit();

    BaseUnit(string name);

    BaseUnit(const UnitStats& baseStats);

    BaseUnit(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max);

    BaseUnit(BaseUnit const& baseUnit);


    pair<double,double> computeFriendForce(BaseUnit const & other);

    pair<double,double> computeEnemyForce(BaseUnit const & other);

	double computeFriendPot(BaseUnit const & other);
	double computeEnemyPot(BaseUnit const & other);

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



    pair<double,double> getPos() const;

    pair<double,double> getMinPos() const;

    pair<double,double> getMaxPos() const;

    double getX() const;

    double getY() const;

    void setPos(double const x, double const y);

    void setPos(const pair<double,double> pos);

    void setMinPos(pair<double,double> const pos);

    void setMaxPos(pair<double,double> const pos);

    void setX(double const x);

    void setY(double const y);

    void setFriendForce(function<pair<double,double>(BaseUnit const &, BaseUnit const &)> func);

    void setEnemyForce(function<pair<double,double>(BaseUnit const &, BaseUnit const &)> func);

    double computeRange(BaseUnit const& other) const;

    pair<double,double> normVecSafe(pair<double,double>const& vec, double const norm) const;

	template <typename T, typename U> double getPot(PlayerState<T>& own, PlayerState<U>& other, std::vector<std::vector<double>> &v)
	{
		double me = 0.0f;
		for(auto buddy :  own.unitList)
		{
		//	me += mFriendFuncPot(*this,*buddy);
			//o << buddy->getPos().first << "\t" << buddy->getPos().second << "\t" << mFriendFuncPot(*this,*buddy)<< std::endl;
		//	std::cout << "\t" << buddy->getPos().first << " "  << buddy->getPos().second << " " << mFriendFuncPot(*this,*buddy) << std::endl; 
			//v.at(buddy->getPos().first).at(buddy->getPos().second) += mFriendFuncPot(*this,*buddy, v);
			me += mFriendFuncPot(*this,*buddy, v, 1);//own.unitList.size());
		//	me += computeFriendPot(*buddy);
		}
		for(auto enemy :  other.unitList)
		{
			//me += mEnemyFuncPot(*this, *enemy);
		//	std::cout << enemy->getPos().first << " " << enemy->getPos().second << " " << mEnemyFuncPot(*this,*enemy) << std::endl;
		//	v.at(enemy->getPos().first).at(enemy->getPos().second) += mEnemyFuncPot(*this,*enemy,v);
			me+=mEnemyFuncPot(*this,*enemy,v, 1);//other.unitList.size());
		//	o << enemy->getPos().first << "\t" << enemy->getPos().second << "\t" << mEnemyFuncPot(*this,*enemy)<< std::endl;
		//	me += computeEnemyPot(*enemy);
		}
		//me = std::min(1e+07, me);
		return std::max(-1e8, std::min(1e8,me));
	}

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
            mPath.push_back(std::make_pair(static_cast<float>(mPos.first),static_cast<float>(mPos.second)));
        }
        if(own.movementTimer <= 0)
        {
            pair<double,double> force(0.0,0.0);
            for(const auto& pot : own.potentialList)
            {
                const pair<double,double> generatedForce = pot.computeForce(*this);
                force.first += generatedForce.first;
                force.second += generatedForce.second;
            }
            for(auto buddy :  own.unitList)
            {
                if(buddy->getHealth() < EPS)
                {
                    continue;
                }

                const pair<double,double> generatedForce = this->computeFriendForce(*buddy);
                force.first += generatedForce.first;
                force.second += generatedForce.second;
            }
            for(auto enemy :  other.unitList)
            {
                if(enemy->getHealth() < EPS)
                {
                    continue;
                }

                const pair<double,double> generatedForce = this->computeEnemyForce(*enemy);
                force.first += generatedForce.first;
                force.second += generatedForce.second;
            }
            double const norm = sqrt(pow(force.first,2)+pow(force.second,2));
            currentForce = this->normVecSafe(force,norm);
        }


        setX(currentForce.first*getSpeed()+getX());
        setY(currentForce.second*getSpeed()+getY());

        // Collision detection, include it again if needed
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
                pair<double,double> shift = normVecSafe(pair<double,double>(x, y), dist);
                setX(getX()+len*shift.first);
                setY(getY()+len*shift.second);
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
                pair<double,double> shift = normVecSafe(pair<double,double>(x, y), dist);
                setX(getX()+len*shift.first);
                setY(getY()+len*shift.second);
            }

        }
        */

    }

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
            if(enemy->getHealth() < EPS)
            {
                continue;
            }
            double const damage = computeDamage(*enemy).total;
            bool const newKill = 0.75 * damage > enemy->getHealth() + enemy->getShield();
            bool const higher = damage > maxDamage;
            if(newKill)
            {
                if(kill)
                {
                    if(higher)
                    {
                        maxDamage = damage;
                        aim = enemy;
                    }
                }
                else
                {
                    maxDamage = damage;
                    aim = enemy;
                }
            }
            else
            {
                if(higher)
                {
                    maxDamage = damage;
                    aim = enemy;
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
    vector<pair<float,float>> getPath() const;
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
};


/* Race Base Unit Types, as there exist unique unit characteristics for the three different races
 f.e. Zerg Units regenerate health and Protoss Units have a regenerating shield */
class TerranUnit : public BaseUnit
{
public:
    TerranUnit();

    TerranUnit(string name);

    TerranUnit(const UnitStats& baseStats);

    TerranUnit(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max);

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

    ZergUnit(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max);

    ZergUnit(ZergUnit const& zergUnit);

    ZergUnit(BaseUnit const& baseUnit);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer >= 1000)
        {
            if(mStats.health > EPS && mStats.health < mStats.maxHealth)
            {
                BaseUnit::addHealth(0.27);
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

    ProtossUnit(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max);

    ProtossUnit(ProtossUnit const& protossUnit);

    ProtossUnit(BaseUnit const& baseUnit);

    void subShield(double const value);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer >= 1000)
        {
            if(mShieldRegenCount == 0 && mStats.shield < mStats.maxShield && mStats.health > EPS)
            {
                BaseUnit::addShield(2.0);
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

    Reaper(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max);

    Reaper(Reaper const& Reaper);

    Reaper(BaseUnit const& baseUnit);

    Reaper(TerranUnit const& terranUnit);

    void subHealth(double const value);

    template <typename T> void regenerate(PlayerState<T>& state)
    {
        if(state.regenerationTimer >= 1000)
        {
            if(mHealthRegenCount == 0 && mStats.health < mStats.maxHealth && mStats.health > EPS)
            {
                BaseUnit::addHealth(2.0);
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


