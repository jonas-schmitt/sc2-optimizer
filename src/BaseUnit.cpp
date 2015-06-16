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


BaseUnit::BaseUnit()
{
    currentForce = Vec2D(0.0);
    for(int i = 0; i < 18; ++i)
    {
        mPossibleDamage[i].valid = false;
    }
}

BaseUnit::BaseUnit(string name) : BaseUnit()
{
    mName = name;
}

BaseUnit::BaseUnit(UnitStats const& baseStats) : BaseUnit()
{
    mStats = baseStats;
}

BaseUnit::BaseUnit(UnitStats const& baseStats, Vec2D min, Vec2D max)
    : BaseUnit(baseStats)
{
    mMinPos = min;
    mMaxPos = max;
}


/* Computes the force that is generated by this BaseUnit */
Vec2D BaseUnit::computeFriendForce(BaseUnit & other)
{
    return mFriendFunc(*this,other);
}

Vec2D BaseUnit::computeEnemyForce(BaseUnit & other)
{
    return mEnemyFunc(*this,other);
}


void BaseUnit::setName(string name)
{
    mName = name;
}

UnitStats& BaseUnit::accessStats()
{
    return this->mStats;
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



void BaseUnit::setMaxHealth (double const value)
{
    mStats.maxHealth = value;
}

void BaseUnit::incMaxHealth (double const value)
{
    mStats.maxHealth += value;
}

void BaseUnit::incHealth (double const value)
{
    double const health_new = mStats.health + value;
    mStats.health = health_new < mStats.maxHealth ? health_new : mStats.maxHealth;
}

double BaseUnit::getMinerals() const
{
    return mStats.minerals;
}

double BaseUnit::getGas() const
{
    return mStats.gas;
}

double BaseUnit::getGroundRange() const
{
    return mStats.groundRange;
}

double BaseUnit::getAirRange() const
{
    return mStats.airRange;
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

double BaseUnit::getSumMaxHealthShield() const
{
    return mStats.sumMaxHealthAndShield;
}

double BaseUnit::getArmor() const
{
    return mStats.armor;
}


double BaseUnit::getSpeed() const
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


bool BaseUnit::isAirUnit() const
{
    return mStats.airUnit;
} //returns true, if its an air unit

double BaseUnit::getSize() const
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
        mStats.health = 0.0;
        return;
    }
    mStats.health = mStats.health > value ? mStats.health - value : 0.0;
}

void BaseUnit::setHealth(double const value)
{
    mStats.health = value;
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

Vec2D BaseUnit::getPos() const
{
    return mPos;
}

Vec2D BaseUnit::getMinPos() const
{
    return mMinPos;
}

Vec2D BaseUnit::getMaxPos() const
{
    return mMaxPos;
}

double BaseUnit::getX() const
{
    return mPos.x;
}

double BaseUnit::getY() const
{
    return mPos.y;
}

void BaseUnit::setPos(double const x, double const y)
{
    setX(x);
    setY(y);
}

void BaseUnit::setPos(const Vec2D pos)
{
    setX(pos.x);
    setY(pos.y);
}

void BaseUnit::setX(double const x)
{
    if(x < mMinPos.x)
    {
        mPos.x = mMinPos.x;
    }
    else if(x > mMaxPos.x)
    {
        mPos.x = mMaxPos.x;
    }
    else
    {
        mPos.x = x;
    }
}

void BaseUnit::setY(double const y)
{

    if(y < mMinPos.y)
    {
        mPos.y = mMinPos.y;
    }
    else if(y > mMaxPos.y)
    {
        mPos.y = mMaxPos.y;
    }
    else
    {
        mPos.y = y;
    }
}


void BaseUnit::setMinPos(Vec2D const pos)
{
    mMinPos = pos;
}

void BaseUnit::setMaxPos(Vec2D const pos)
{
    mMaxPos = pos;
}

void BaseUnit::setFriendForce(function<Vec2D(BaseUnit &, BaseUnit &)> func)
{
    mFriendFunc = func;
}
void BaseUnit::setEnemyForce(function<Vec2D(BaseUnit &, BaseUnit &)> func)
{
    mEnemyFunc = func;
}


double BaseUnit::getResources() const
{
    return mStats.minerals+GASTOMINERALS*mStats.gas;
}


double BaseUnit::getAirAttack() const
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

double BaseUnit::getGroundAttack() const
{
    return mStats.groundAttack;
}

double BaseUnit::getGAUpgradeBonus() const
{
    return mStats.gaUpgrade;
}
double BaseUnit::getAAUpgradeBonus() const
{
    return mStats.aaUpgrade;
}

double BaseUnit::getArmorUpgradeBonus() const
{
    return mStats.armorUpgradeBonus;
}

int BaseUnit::getAttackUpgrade() const
{
    return mAttackUpgrade;
}

int BaseUnit::getArmorUpgrade() const
{
    return mArmorUpgrade;
}

vector<Attribute> const& BaseUnit::getAttributes() const
{
    return mStats.attributes;
}

vector<Bonus> const& BaseUnit::getBonuses() const
{
    return mStats.bonuses;
}


int BaseUnit::getShieldUpgrade() const
{
    return mShieldUpgrade;
}

void BaseUnit::setShieldUpgrade(int value)
{
    mShieldUpgrade = value;
}

double BaseUnit::getAttackMultiplier() const
{
    return mAttackMultiplier;
}

void BaseUnit::setAttackMultiplier(double value)
{
    mAttackMultiplier = value;
}

double BaseUnit::getDefenseMultiplier() const
{
    return mDefenseMultiplier;
}

void BaseUnit::setDefenseMultiplier(double value)
{
    mDefenseMultiplier = value;
}

double BaseUnit::getDefenseSubtractor() const
{
    return mDefenseSubtractor;
}

void BaseUnit::setDefenseSubtractor(double value)
{
    mDefenseSubtractor = value;
}



int BaseUnit::getIdentifier() const
{
    return mId;
}

void BaseUnit::setIdentifier(int value)
{
    mId = value;
}

Vec2D BaseUnit::getStartPos() const
{
    return mStartPos;
}

void BaseUnit::setStartPos(const Vec2D &value)
{
    mStartPos = value;
}

void BaseUnit::resetPos()
{
    mPos = mStartPos;
}


int BaseUnit::getMovementTimer() const
{
    return mMovementTimer;
}

void BaseUnit::setMovementTimer(int value)
{
    mMovementTimer = value;
}

int BaseUnit::getMovementUpdate() const
{
    return mMovementUpdate;
}

void BaseUnit::setMovementUpdate(int value)
{
    mMovementUpdate = value;
}

void BaseUnit::resetTimer ()
{
    mMovementTimer = 0;
    mAttackTimer = 0;
}

void BaseUnit::decMovementTimer()
{
    if(mMovementTimer > 0)
    {
        mMovementTimer -= mTimeSlice;
    }
}


double BaseUnit::getMoveDist() const
{
    return mMoveDist;
}

void BaseUnit::multSpeed (double value)
{
    mStats.speed *= value;
    mMoveDist *= value;
}

bool BaseUnit::isDead () const
{
    return mStats.health < EPS;
}


size_t BaseUnit::getChromosomeStartPosition() const
{
    return mChromosomeStartPosition;
}
Damage BaseUnit::computeDamage(BaseUnit const& other) const
{
    double totalDamage = 0.0;
    if(other.isAirUnit())
    {
        totalDamage += mAttackUpgrade * mStats.aaUpgrade + mStats.airAttack;
    }
    else
    {
        totalDamage += mAttackUpgrade * mStats.gaUpgrade + mStats.groundAttack;
    }
    for(Bonus const& bonus : mStats.bonuses)
    {
        vector<Attribute> const& attributes = other.getAttributes();
        if(std::includes(attributes.begin(), attributes.end(), bonus.attributes.begin(), bonus.attributes.end()))
        {
            totalDamage += mAttackUpgrade * bonus.upgrade + bonus.base;
        }
    }
    double const shield = other.getShield();
    totalDamage *= mAttackMultiplier;
    if(shield > 0)
    {
        totalDamage -= other.getShieldUpgrade() * other.getDefenseMultiplier();
    }
    totalDamage -= other.getDefenseSubtractor();
    double shieldDamage;
    double healthDamage = totalDamage - shield;
    if(healthDamage > 0)
    {
        shieldDamage = shield;
        healthDamage -= (other.getArmorUpgrade() * other.getArmorUpgradeBonus() + other.getArmor()) * other.getDefenseMultiplier();
        if(healthDamage < 0.0)
        {
            healthDamage = 0.0;
        }
        totalDamage = shieldDamage + healthDamage;
    }
    else
    {
        shieldDamage = totalDamage;
        healthDamage = 0.0;
    }

    Damage damage(shieldDamage, healthDamage, totalDamage);

    if(damage.total < 0.5)
    {
        damage.total = 0.5;
        if(shield > 0.5)
        {
            damage.shield = 0.5;
            damage.health = 0.0;
        }
        else
        {
            damage.shield = shield;
            damage.health = 0.5 - shield;
        }
    }

    return damage;

}

bool BaseUnit::attackPossible(BaseUnit const& enemy)
{
    if(enemy.getHealth() < EPS)
    {
        return false;
    }
    double const x = enemy.getX() - this->getX();
    double const y = enemy.getY() - this->getY();
    double dist = std::sqrt(x*x+y*y);
    if(std::isnan(dist))
    {
        dist = 0.0;
    }
    return computeRange (enemy) > dist;

}

bool BaseUnit::attack(BaseUnit& unit)
{
    Damage const& damage = mPossibleDamage[unit.getIdentifier()];
    if(damage.shield > 0)
    {
        unit.subShield(damage.shield);
    }
    if(damage.health > 0)
    {
        unit.subHealth(damage.health);
    }
    if(unit.isAirUnit())
    {
        this->setAttackTimer(this->getAACooldown());
    }
    else
    {
        this->setAttackTimer(this->getGACooldown());
    }
    if(unit.getHealth() < EPS)
    {
        return true;
    }
    else
    {
        mTarget = &unit;
        return false;
    }
}


double BaseUnit::getMaxDist() const
{
    double const x = this->mMaxPos.x-this->mMinPos.x;
    double const y = this->mMaxPos.y-this->mMinPos.y;
    return std::sqrt(x*x+y*y);
}


double BaseUnit::computeRange(BaseUnit const& other) const
{
    return other.isAirUnit() ? computeAirRange (other) : computeGroundRange (other);
}

//TODO These functions must return 0 if the unit can generally not attack the respective type
double BaseUnit::computeGroundRange(BaseUnit const& other) const
{
    return getGroundRange() + getSize() + other.getSize();
}
double BaseUnit::computeAirRange(BaseUnit const& other) const
{
    return getAirRange() + getSize() + other.getSize();
}

double BaseUnit::getPhenotype(size_t const pos) const
{
    return mPhenotype[pos];
}

void BaseUnit::setChromosomeStartPosition(size_t const pos)
{
    mChromosomeStartPosition = pos;
}

void BaseUnit::setChromosome(Chromosome const& chromosome)
{

    computeTemporaryValues();
    mChromosome = chromosome.data();
    mPhenotype.reserve(mNGenes);
    double constexpr max = 1.0/(std::pow(2, NBITS) - 1);
    for(int i = 0; i < chromosome.size(); ++i)
    {
        mPhenotype.push_back(static_cast<double>(mChromosome[mChromosomeStartPosition + i].to_ulong())*max);
    }
}

void BaseUnit::setChromosome(Chromosome const& chromosome, size_t const startPos)
{
    mChromosomeStartPosition = startPos;
    setChromosome(chromosome);
}

void BaseUnit::setTracking(bool const tracking)
{
    mTracking = tracking;
}

void BaseUnit::reservePathStorage(size_t const sz)
{
    mPath.reserve(sz);
}

vector<Vec2Df> BaseUnit::getPath() const
{
    return mPath;
}

void BaseUnit::clearPath()
{
    mPath.clear();
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
    if(mAttackTimer > 0)
    {
        mAttackTimer -= mTimeSlice;
    }

}


int BaseUnit::getTimeSlice() const
{
    return mTimeSlice;
}

void BaseUnit::setTimeSlice(int value)
{
    mTimeSlice = value;
}

void BaseUnit::computeTemporaryValues()
{
    mMoveDist = getSpeed()*mTimeSlice*1e-3;
    mMovementUpdateBackup = mMovementUpdate;
}

void BaseUnit::initUpgrades (vector<int> const& flags)
{
    mAttackUpgrade = flags[0];
    mArmorUpgrade = flags[1];
}



