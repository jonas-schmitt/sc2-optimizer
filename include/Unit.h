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
using std::function;
using std::pow;
using std::sqrt;
using std::list;
using std::string;
using std::vector;

enum class AttackModus
{
    splash, projectile, linear
};

enum class Attribute
{
    light, armored, biological, psyonic, massive, structure
};


struct UnitStats
{
    UnitStats()
    {
        minerals = 0;
        gas = 0;
        gdps = 0;
        adps = 0;
        groundRange = 0;
        airRange = 0;
        health = 0;
        maxHealth = 0;
        shield = 0;
        maxShield = 0;
        armor = 0;
        shieldArmor = 0;
        sight = 0;
        speed = 0;
        acceleration = 0;
        energy = 0;
        maxEnergy = 0;
        size = 0;
        airUnit = false;
    }

    UnitStats(const UnitStats& stats)
        : minerals(stats.minerals), gas(stats.gas), gdps(stats.gdps), adps(stats.adps), groundRange(stats.groundRange),
          airRange(stats.airRange), health(stats.health), maxHealth(stats.maxHealth), shield(stats.shield), maxShield(stats.maxShield), armor(stats.armor), shieldArmor(stats.shieldArmor), sight(stats.sight),
          speed(stats.speed), acceleration(stats.acceleration), energy(stats.energy), maxEnergy(stats.maxEnergy), size(stats.size), airUnit(stats.airUnit), attributes(stats.attributes),
          bonusDps(stats.bonusDps), attackMods(stats.attackMods)
    {};

    float minerals;
    float gas;
    float gdps;
    float adps;
    float groundRange;
    float airRange;
    double health;
    double maxHealth;
    double shield;
    double maxShield;
    float armor;
    float shieldArmor = 0;
    float sight;
    float speed;
    float acceleration;
    double energy;
    double maxEnergy;
    float size;
    bool airUnit;
    vector<Attribute> attributes;
    vector<pair<float,Attribute> > bonusDps;
    vector<AttackModus> attackMods; 
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
        if(dist > ownRange)
        {
			val = own.getXGene(1)*enemy.getResources()+own.getXGene(2)*own.computeDps(enemy)+own.getYGene(2);
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
            val2 = (-own.getYGene(6)-own.getXGene(3)*own.getResources()-own.getXGene(4)*enemy.computeDps(own));
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
        if(dist > ownRange)
        {
            res1 = own.normVecSafe(pair<double,double>(x,y), dist);
            double const val = own.getXGene(1)*enemy.getResources()+own.getXGene(2)*own.computeDps(enemy)+own.getYGene(2);
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
            double const val = -own.getYGene(6)-own.getXGene(3)*own.getResources()-own.getXGene(4)*enemy.computeDps(own);
            res2.first *= val;
            res2.second *= val;
        }
        else
        {
            res2 = pair<double,double>(0.,0.);
        }

        return pair<double,double>(res1.first+res2.first, res1.second+res2.second);
    };

    double computeDamageDealt(BaseUnit const& unit);

    bool attack(BaseUnit* unit);

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

    float getShieldArmor() const;

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

    void incShieldArmor();

    void decShieldArmor();

    void resetHealth();

    void resetShield();

    void resetEnergy();

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

    double computeDps(BaseUnit const& other) const;

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
        pair<double,double> force(0,0);
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
            /*
            double const x = buddy->getX()-this->getX();
            double const y = buddy->getY()-this->getY();
            double  dist = std::sqrt(pow(x,2)+pow(y,2));
            if(this->getHash() != buddy->getHash() && dist < buddy->getSize() + this->getSize() + this->getSpeed())
            {
                mInRange.emplace_back(buddy);

            }*/
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
            /*
            double const x = enemy->getX()-this->getX();
            double const y = enemy->getY()-this->getY();
            double  dist = std::sqrt(pow(x,2)+pow(y,2));
            if(dist < enemy->getSize() + this->getSize() + this->getSpeed())
            {
                mInRange.emplace_back(enemy);

            }*/
            const pair<double,double> generatedForce = this->computeEnemyForce(*enemy);
            force.first += generatedForce.first;
            force.second += generatedForce.second;
        }

        currentForce.first = fmin(fmax(currentForce.first+force.first,-LIMIT),LIMIT);
        currentForce.second = fmin(fmax(currentForce.second+force.second,-LIMIT),LIMIT);
        double norm = sqrt(pow(currentForce.first,2)+pow(currentForce.second,2));

        auto res = this->normVecSafe(currentForce,norm);
        setX(res.first*getSpeed()+getX());
        setY(res.second*getSpeed()+getY());
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

    }

    template<typename T> void attack(PlayerState<T>& other)
    {
        if(this->getHealth() < EPS)
        {
            return;
        }
        if(other.unitList.empty())
        {
            return;
        }

        auto aim = *(other.unitList.begin());
        double maxDmg = computeDamageDealt(*aim);
        double dmg = 0;
        double wastedDmg = 0;
        bool kill = false;
        for(auto enemy : other.unitList)
        {
            if(enemy->getHealth() < EPS)
            {
                continue;
            }

            dmg = computeDamageDealt(*enemy);
            if(3./4.*dmg > enemy->getHealth()+enemy->getShield())
            {
                if(kill)
                {
                    if(dmg-enemy->getHealth()-enemy->getShield() < wastedDmg)
                    {
                        wastedDmg = dmg-enemy->getHealth()-enemy->getShield();
                        maxDmg = dmg;
                        aim = enemy;
                    }
                }
                else
                {
                    wastedDmg = dmg-enemy->getHealth();
                    maxDmg = dmg;
                    kill = true;
                    aim = enemy;
                }
            }
            else if(!kill && dmg > maxDmg)
            {
                maxDmg = dmg;
                aim = enemy;
            }
        }
        if(maxDmg < EPS)
        {
            return;
        }
        if(attack(aim))
        {
            --other.unitCount;
        }
    }

    int getXGene(size_t pos) const;
    int getYGene(size_t pos) const;
    void setGenes(UnitGenes const& genes);
    size_t getHash() const;

    void setTracking(bool const tracking);
    void reservePathStorage(size_t const sz);
    vector<pair<float,float>> getPath() const;
    void clearPath();

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

    void regenerate();
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

    void regenerate();

};

class ProtossUnit : public BaseUnit
{
protected:
    int shieldCount;
public:
    ProtossUnit();

    ProtossUnit(string name);

    ProtossUnit(const UnitStats& baseStats);

    ProtossUnit(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max);

    ProtossUnit(ProtossUnit const& protossUnit);

    ProtossUnit(BaseUnit const& baseUnit);

    void subShield(double const value);

    void regenerate();
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
    int regenCount;
public:
    Reaper();

    Reaper(string name);

    Reaper(const UnitStats& baseStats);

    Reaper(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max);

    Reaper(Reaper const& Reaper);

    Reaper(BaseUnit const& baseUnit);

    Reaper(TerranUnit const& terranUnit);

    void subHealth(double const value);

    void regenerate();

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


