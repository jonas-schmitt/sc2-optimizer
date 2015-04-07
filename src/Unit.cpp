#include<utility>
#include<functional>
#include<cmath>
#include<list>

#include "../include/Unit.h"

using std::pair;
using std::function;
using std::pow;
using std::sqrt;
using std::list;
using std::string;


BaseUnit::BaseUnit(){}

BaseUnit::BaseUnit(string name) : mName(name)
{
	currentForce = std::make_pair<double,double>(0,0);
}

BaseUnit::BaseUnit(UnitStats const& baseStats) : mStats(baseStats)
{
	currentForce = std::make_pair<double,double>(0,0);
}

BaseUnit::BaseUnit(UnitStats const& baseStats, pair<double,double>min, pair<double,double>max)
    : mStats(baseStats), mMinPos(min), mMaxPos(max)
{
	currentForce = std::make_pair<double,double>(0,0);
}

BaseUnit::BaseUnit(BaseUnit const& baseUnit) : mStats(baseUnit.mStats), mPos(baseUnit.mPos), mMinPos(baseUnit.mMinPos), mMaxPos(baseUnit.mMaxPos),
    mFriendFunc(baseUnit.mFriendFunc), mEnemyFunc(baseUnit.mEnemyFunc)
{
	currentForce = std::make_pair<double,double>(0,0);
}

/* Computes the force that is generated by this BaseUnit */
pair<double,double> BaseUnit::computeFriendForce(BaseUnit const& other)
{
    return mFriendFunc(*this,other);
}

pair<double,double> BaseUnit::computeEnemyForce(BaseUnit const& other)
{
    return mEnemyFunc(*this,other);
}

double BaseUnit::computeFriendPot(BaseUnit const& other)
{
	//return mFriendFuncPot(*this, other);
	return 0;
}

double BaseUnit::computeEnemyPot(BaseUnit const& other)
{
	//return mEnemyFuncPot(*this, other);
	return 0;
}

void BaseUnit::setName(string name)
{
    mName = name;
}

string BaseUnit::getName() const
{
    return mName;
}

// getter and setter
void BaseUnit::setStats(const UnitStats& newStats)
{
    this->mStats = newStats;
}

float BaseUnit::getMinerals() const
{
    return mStats.minerals;
}

float BaseUnit::getGas() const
{
    return mStats.gas;
}

float BaseUnit::getGdps() const
{
    return mStats.gdps;
}

float BaseUnit::getAdps() const
{
    return mStats.adps;
}

float BaseUnit::getGroundRange() const
{
    return mStats.groundRange+mStats.size;
}

float BaseUnit::getAirRange() const
{
    return mStats.airRange+mStats.size;
}

double BaseUnit::getHealth() const
{
    return mStats.health;
}

double BaseUnit::getMaxHealth() const
{
    return mStats.maxHealth;
}

double BaseUnit::getShield() const
{
    return mStats.shield;
}

double BaseUnit::getMaxShield() const
{
    return mStats.maxShield;
}

float BaseUnit::getArmor() const
{
    return mStats.armor;
}


float BaseUnit::getSpeed() const
{
    return mStats.speed;
}


double BaseUnit::getEnergy() const
{
    return mStats.energy;
}

double BaseUnit::getMaxEnergy() const
{
    return mStats.maxEnergy;
}

size_t BaseUnit::getHash() const
{
    return mGenes.getHash();
}

bool BaseUnit::isAirUnit() const
{
    return mStats.airUnit;
} //returns true, if its an air unit

float BaseUnit::getSize() const
{
	return mStats.size;
}


void BaseUnit::addHealth(double const value)
{
    if(mStats.health + value > mStats.maxHealth)
    {
        mStats.health = mStats.maxHealth;
    }
    else
    {
        mStats.health += value;
    }
}

void BaseUnit::subHealth(double const value)
{
    if(mStats.health < EPS)
    {
        mStats.health = 0;
        return;
    }
    double newHealth = mStats.health - value;
    if(newHealth < EPS)
    {
        mStats.health = 0.;
    }
    else
    {
        mStats.health = newHealth;
    }
}

void BaseUnit::addShield(double const value)
{
    if(mStats.shield + value > mStats.maxShield)
    {
        mStats.shield = mStats.maxShield;
    }
    else
    {
        mStats.shield += value;
    }
}

void BaseUnit::subShield(double const value)
{
    if(mStats.shield < EPS)
    {
        mStats.shield = 0;
        return;
    }
    mStats.shield = mStats.shield > value ? mStats.shield - value : 0;
}

void BaseUnit::resetHealth()
{
    mStats.health = mStats.maxHealth;
}

void BaseUnit::resetShield()
{
    mStats.shield = mStats.maxShield;
}

void BaseUnit::resetEnergy()
{
    mStats.energy = mStats.maxEnergy;
}

void BaseUnit::incArmor()
{
    ++mStats.armor;
}

void BaseUnit::decArmor()
{
    --mStats.armor;
}

pair<double,double> BaseUnit::getPos() const
{
    return mPos;
}

pair<double,double> BaseUnit::getMinPos() const
{
    return mMinPos;
}

pair<double,double> BaseUnit::getMaxPos() const
{
    return mMaxPos;
}

double BaseUnit::getX() const
{
    return mPos.first;
}

double BaseUnit::getY() const
{
    return mPos.second;
}

void BaseUnit::setPos(double const x, double const y)
{
    setX(x);
    setY(y);
}

void BaseUnit::setPos(const pair<double,double> pos)
{
    setX(pos.first);
    setY(pos.second);
}

void BaseUnit::setX(double const x)
{
    if(x < mMinPos.first)
    {
        mPos.first = mMinPos.first;
    }
    else if(x > mMaxPos.first)
    {
        mPos.first = mMaxPos.first;
    }
    else
    {
        mPos.first = x;
    }
}

void BaseUnit::setY(double const y)
{

    if(y < mMinPos.second)
    {
        mPos.second = mMinPos.second;
    }
    else if(y > mMaxPos.second)
    {
        mPos.second = mMaxPos.second;
    }
    else
    {
        mPos.second = y;
    }
}


void BaseUnit::setMinPos(pair<double,double> const pos)
{
    mMinPos = pos;
}

void BaseUnit::setMaxPos(pair<double,double> const pos)
{
    mMaxPos = pos;
}

void BaseUnit::setFriendForce(function<pair<double,double>(BaseUnit const&, BaseUnit const&)> func)
{
    mFriendFunc = func;
}
void BaseUnit::setEnemyForce(function<pair<double,double>(BaseUnit const&, BaseUnit const&)> func)
{
    mEnemyFunc = func;
}


float BaseUnit::getResources() const
{
    return mStats.minerals+GASTOMINERALS*mStats.gas;
}


float BaseUnit::getAirAttack() const
{
    return mStats.airAttack;
}


int BaseUnit::getGACooldown() const
{
    return mStats.gaCooldown;
}


int BaseUnit::getAACooldown() const
{
    return mStats.aaCooldown;
}

float BaseUnit::getGroundAttack() const
{
    return mStats.groundAttack;
}

float BaseUnit::getGAUpgrade() const
{
    return mStats.gaUpgrade;
}
float BaseUnit::getAAUpgrade() const
{
    return mStats.aaUpgrade;
}

float BaseUnit::getArmorUpgrade() const
{
    return mStats.armorUpgrade;
}


double BaseUnit::computeDamageDealt(BaseUnit const& unit)
{
    if(unit.getHealth() < EPS)
    {
        return 0;
    }
    double&& x = unit.getX() - this->getX();
    double&& y = unit.getY() - this->getY();
    double dist = std::sqrt(pow(x,2)+pow(y,2));
    if(std::isnan(dist))
    {
        dist = 0.;
    }
    double result = 0;
    if(unit.isAirUnit() && this->getAirRange() > dist)
    {
        result = this->getAirAttack()-unit.getArmor();
    }
    else if(!unit.isAirUnit() && this->getGroundRange() > dist)
    {
        result = this->getGroundAttack()-unit.getArmor();
    }
    return result < EPS ? 0 : result;
}

bool BaseUnit::attack(BaseUnit* unit)
{
    if(this->getHealth() < EPS || unit->getHealth() < EPS)
    {
        return false;
    }

    double&& x = unit->getX() - this->getX();
    double&& y = unit->getY() - this->getY();
    double dist = std::sqrt(pow(x,2)+pow(y,2));
    if(std::isnan(dist))
    {
        dist = 0.;
    }
    if((unit->isAirUnit() && dist > this->getAirRange()) || (!unit->isAirUnit() && dist > this->getGroundRange()))
    {
        return false;
    }
    double damage = unit->isAirUnit() ? this->getAdps() : this->getGdps();
    if(damage < EPS)
    {
        return false;
    }
    damage -= unit->getArmor();
    if(damage-0.5 < EPS)
    {
        damage = 0.5;
    }
    double oldShield = unit->getShield();
    unit->subShield(damage);
    if(damage-oldShield < EPS)
    {
        return false;
    }
    damage -= unit->getArmor();
    if(damage-0.5 < EPS)
    {
        damage = 0.5;
    }
    bool kill = false;
    if(unit->getHealth() > EPS && damage > unit->getHealth())
    {
        kill = true;
    }
    if(unit->isAirUnit())
    {
        this->setAttackTimer(this->getAACooldown());
    }
    else
    {
        this->setAttackTimer(this->getGACooldown());
    }
    unit->subHealth(damage);
    return kill;
}

pair<double,double> BaseUnit::normVecSafe(pair<double,double>const& vec, double const norm) const
{
    pair<double,double> result;
    if(std::isinf (vec.first/norm) || std::isnan(vec.first/norm) || std::isinf (vec.second/norm) || std::isnan(vec.second/norm) || norm < EPS)
    {
        // choose direction that brings the unit as much away from the border as possible

        if(vec.first > EPS)
        {
            result.first = STDLEN;
        }
        else if(vec.first < -EPS)
        {
            result.first = -STDLEN;
        }
        else
        {
            result.first = 0;
        }
        if(vec.second > EPS)
        {
            result.second = STDLEN;
        }
        else if(vec.second < -EPS)
        {
            result.second = -STDLEN;
        }
        else
        {
            result.second = 0;
        }
    }
    else
    {
        result.first = vec.first/norm;
        result.second = vec.second/norm;
    }
    return result;
}

double BaseUnit::getMaxDist() const
{
    double const x = this->mMaxPos.first-this->mMinPos.first;
    double const y = this->mMaxPos.second-this->mMinPos.second;
    return std::sqrt(pow(x,2)+pow(y,2));
}
double BaseUnit::computeDamage(BaseUnit const& other) const
{
    return other.isAirUnit() ? this->getAdps()-other.getArmor() : this->getGdps()-other.getArmor();
}

double BaseUnit::computeRange(BaseUnit const& other) const
{
    return other.isAirUnit() ? this->getAirRange() : this->getGroundRange();
}

int BaseUnit::getXGene(size_t pos) const
{
    return mGenes.getX(pos);
}
int BaseUnit::getYGene(size_t pos) const
{
    return mGenes.getY(pos);
}

void BaseUnit::setGenes(UnitGenes const& genes)
{
    mGenes = genes;
}

void BaseUnit::setTracking(bool const tracking)
{
    mTracking = tracking;
}

void BaseUnit::reservePathStorage(size_t const sz)
{
    mPath.reserve(sz);
}

vector<pair<float,float>> BaseUnit::getPath() const
{
    return mPath;
}

void BaseUnit::clearPath()
{
    mPath.clear();
}

int BaseUnit::getMovementTimer() const
{
    return mMovementTimer;
}

void BaseUnit::setMovementTimer(int value)
{
    mMovementTimer = value;
}

void BaseUnit::decMovementTimer()
{
    if(mMovementTimer < mTimeSlice)
    {
        mMovementTimer = 0;
        return;
    }
    mMovementTimer -= mTimeSlice;
}

void BaseUnit::resetMovementTimer()
{
    mMovementTimer = mMovementUpdate;
}

int BaseUnit::getMovementUpdate() const
{
    return mMovementUpdate;
}

void BaseUnit::setMovementUpdate(int value)
{
    mMovementUpdate = value;
}


int BaseUnit::getAttackTimer() const
{
    return mAttackTimer;
}
void BaseUnit::setAttackTimer(int value)
{
    mAttackTimer = value;
}

void BaseUnit::decAttackTimer()
{
    if(mAttackTimer < mTimeSlice)
    {
        mAttackTimer = 0;
        return;
    }
    mAttackTimer -= mTimeSlice;
}


int BaseUnit::getTimeSlice() const
{
    return mTimeSlice;
}

void BaseUnit::setTimeSlice(int value)
{
    mTimeSlice = value;
}




TerranUnit::TerranUnit()
    : BaseUnit()
{}

TerranUnit::TerranUnit(string name)
    : BaseUnit(name)
{}

TerranUnit::TerranUnit(const UnitStats& baseStats)
    : BaseUnit(baseStats)
{}

TerranUnit::TerranUnit(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max)
    : BaseUnit(baseStats,min,max)
{}

TerranUnit::TerranUnit(TerranUnit const& terranUnit)
 : BaseUnit(terranUnit)
{}

TerranUnit::TerranUnit(BaseUnit const& baseUnit)
    : BaseUnit(baseUnit)
{}

void TerranUnit::regenerate()
{}




ZergUnit::ZergUnit()
    : BaseUnit()
{}

ZergUnit::ZergUnit(string name)
    : BaseUnit(name)
{}

ZergUnit::ZergUnit(const UnitStats& baseStats)
    : BaseUnit(baseStats)
{}

ZergUnit::ZergUnit(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max)
    : BaseUnit(baseStats,min,max)
{}

ZergUnit::ZergUnit(ZergUnit const& zergUnit)
 : BaseUnit(zergUnit)
{}

ZergUnit::ZergUnit(BaseUnit const& baseUnit)
    : BaseUnit(baseUnit)
{}

void ZergUnit::regenerate()
{
    if(mStats.health > EPS && mStats.health < mStats.maxHealth)
    {
        BaseUnit::addHealth(0.27);
    }
}









ProtossUnit::ProtossUnit()
    : BaseUnit(), shieldCount(0)
{}

ProtossUnit::ProtossUnit(string name)
    : BaseUnit(name), shieldCount(0)
{}

ProtossUnit::ProtossUnit(const UnitStats& baseStats)
    : BaseUnit(baseStats), shieldCount(0)
{}

ProtossUnit::ProtossUnit(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max)
    : BaseUnit(baseStats,min,max), shieldCount(0)
{}

ProtossUnit::ProtossUnit(ProtossUnit const& protossUnit)
 : BaseUnit(protossUnit), shieldCount(0)
{}

ProtossUnit::ProtossUnit(BaseUnit const& baseUnit)
    : BaseUnit(baseUnit), shieldCount(0)
{}


void ProtossUnit::subShield(double const value)
{
    shieldCount = 10;
    BaseUnit::subShield(value);
}

void ProtossUnit::regenerate()
{
    if(shieldCount == 0 && mStats.shield < mStats.maxShield && mStats.health > EPS)
    {
        BaseUnit::addShield(2.0);
    }
    else
    {
        --shieldCount;
    }

}




Reaper::Reaper()
    : TerranUnit(), regenCount(0)
{}

Reaper::Reaper(string name)
    : TerranUnit(name), regenCount(0)
{}

Reaper::Reaper(const UnitStats& baseStats)
    : TerranUnit(baseStats), regenCount(0)
{}

Reaper::Reaper(const UnitStats& baseStats, pair<double,double>min, pair<double,double>max)
    : TerranUnit(baseStats,min,max), regenCount(0)
{}

Reaper::Reaper(Reaper const& Reaper)
 : TerranUnit(Reaper), regenCount(0)
{}

Reaper::Reaper(BaseUnit const& baseUnit)
    : TerranUnit(baseUnit), regenCount(0)
{}

Reaper::Reaper(TerranUnit const& terranUnit)
    : TerranUnit(terranUnit), regenCount(0)
{}

void Reaper::subHealth(double const value)
{
    regenCount = 10;
    BaseUnit::subHealth(value);
}

void Reaper::regenerate()
{
    if(regenCount == 0 && mStats.health < mStats.maxHealth && mStats.health > EPS)
    {
        BaseUnit::addHealth(2.0);
    }
    else
    {
        --regenCount;
    }

}


