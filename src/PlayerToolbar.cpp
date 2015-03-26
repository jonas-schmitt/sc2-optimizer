#include "PlayerToolbar.h"

PlayerToolbar::PlayerToolbar(QWidget *parent) : QTreeWidget(parent)
{
	setMaximumWidth(200);
	lastSize = -1;
}

PlayerToolbar::~PlayerToolbar()
{
}

void PlayerToolbar::setHeader(const QString name)
{
	setHeaderLabel(name);
}

void PlayerToolbar::setUnitNames(std::vector<std::string> name)
{
	mNames = name;
}

void PlayerToolbar::setUnitStats(std::vector<std::string> stat)
{
	mStats = stat;
}

void PlayerToolbar::setStatValues(std::vector<std::string> value)
{
	mValues = value;
}

void PlayerToolbar::showUnits()
{
	unsigned int numOfStats = mStats.size();
	unsigned int numOfUnits = mNames.size();

	if (lastSize != numOfUnits)
	{
		for (auto in : mStatBox)
			delete in;
		for (auto in : mUnitBox)
		{
			delete in;
		}
		mUnitBox.clear();
		mStatBox.clear();
		mUnitCheck.clear();
	}
	lastSize = numOfUnits;

	if (mUnitBox.isEmpty() == true)
	{
		for (unsigned int i = 0; i < numOfUnits; ++i)
		{
			mUnitBox.append(new QTreeWidgetItem());
			mUnitCheck.append(false);
		}
	}
	if (mStatBox.isEmpty() == true)
	{
		for (unsigned int i = 0; i < numOfUnits; ++i)
		{
			for (unsigned int j = 0; j < numOfStats; ++j)
					mStatBox.append(new QTreeWidgetItem());
		}
	}

	//for number of units in unitlist
	for (unsigned int i = 0; i < numOfUnits; ++i)
	{
		QString name = QString::fromStdString(mNames.at(i));
		mUnitBox[i]->setText(0, name);
		addTopLevelItem(mUnitBox[i]);
		mUnitCheck[i] = false;

		for (unsigned int j = 0; j < numOfStats; ++j)
		{
			QString statName = QString::fromStdString(mStats.at(j));
			statName.append(": ");
			statName.append(QString::fromStdString(mValues.at((i*numOfStats)+j)));
			mStatBox[(i*numOfStats)+j]->setText(0, statName);
			mUnitBox[i]->addChild(mStatBox[(i*numOfStats)+j]);

		}
	}
}

void PlayerToolbar::resetButtonPressed()
{
	//reset colors and marks
	QBrush br;
	br.setColor(Qt::black);
	for (unsigned int i = 0; i < mNames.size(); ++i)
	{
		mUnitBox[i]->setForeground(0, br);
	}
}

void PlayerToolbar::unitMarked(QTreeWidgetItem* item)
{
	if(!item->parent())
	{
		int idx = indexOfTopLevelItem(item);
		//check if unit is already marked
		if (mUnitCheck[idx] == false)
		{
			//unit is not marked, so mark it then
			QBrush br;
			br.setColor(Qt::red);
			mUnitBox[idx]->setForeground(0, br);
			emit(markedUnit(idx));
		}
		else
		{
			//unit should be unselected
			QBrush br;
			br.setColor(Qt::black);
			mUnitBox[idx]->setForeground(0, br);
			emit(markedUnit(idx));
		}
	}
	else
	{
		//don't know what to do here
/*		int idx_p = indexOfTopLevelItem(item->parent());
		int idx_c = item->parent()->indexOfChild(item);
//		std::cout << "Parent: " << idx_p << " Child: " << idx_c << std::endl;
*/
	}
}


