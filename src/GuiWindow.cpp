#include "GuiWindow.h"

static size_t steps=0;

void GuiWindow::threadWrapper()
{
	while (true)
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(150));
		if (updateIntervals <= 0)
		{
			continue;
		}
		if (startFunc(updateIntervals) == true)
		{
			emit showUnits();
			std::cout << "Sim finished" << std::endl;
			return;
		}
		std::cout << steps << std::endl;
		steps += updateIntervals;
		emit showUnits();
	}
}

GuiWindow::GuiWindow(QWidget* parent) : QMainWindow(parent)
{
	mWidget = new QWidget();
	mHLayout = new QHBoxLayout();
	mCentralLayout = new QVBoxLayout();
	mToolbarLayout = new QHBoxLayout();

	mPlayGround = new PlayGround();
	mPlayer1 = new PlayerToolbar();
	mPlayer2 = new PlayerToolbar();

	initGuiElements();

	mWidget->setLayout(mHLayout);
	setCentralWidget(mWidget);
	setWindowTitle("Simulation");

	std::pair<int,int> initPos;
	initPos = std::make_pair(-1,-1);
	mPlayer1BaseCoords = initPos;
	mPlayer2BaseCoords = initPos;
	playable = true;
	updateIntervals = -1;

	makeConnections();
	globalThread = std::thread(&GuiWindow::threadWrapper, this);
	qRegisterMetaType<QVector<int> >("QVector<int>");
}

GuiWindow::~GuiWindow()
{
	if(mPlayGround) delete(mPlayGround);
	if(mPlayer1) 	delete(mPlayer1);
	if(mPlayer2) 	delete(mPlayer2);
}

std::pair<int,int> GuiWindow::getBaseCoordinates(bool player)
{
	if (player == true)
		return mPlayer1BaseCoords;
	else
		return mPlayer2BaseCoords;
}

void GuiWindow::scaleFactor(std::pair<double,double>size)
{
	std::pair<double,double> plSize = mPlayGround->getPlayGroundSize();
	double ratioX = plSize.first/size.first;
	double ratioY = plSize.second/size.second;
	scale = std::make_pair(ratioX,ratioY);
}

template<typename T, typename U>
void GuiWindow::setUnits(T *pl1, U *pl2)
{
    scaleFactor(std::make_pair(pl1->maxPos.x-pl1->minPos.x, pl1->maxPos.y-pl1->minPos.y));
	std::vector<std::pair<double,double>> coord1;
	std::vector<std::pair<double,double>> coord2;
	std::vector<double> sizes;
	std::vector<double> ranges;
	for (auto in : pl1->unitList)
	{
        if(in->getHealth() < EPS)
        {
            continue;
        }
		std::pair<double,double> tmp;
        tmp = std::make_pair(in->getPos().x*scale.first, in->getPos().y*scale.second);
		coord1.push_back(tmp);
		sizes.push_back((in->getSize())*scale.first);
        ranges.push_back((std::max(in->getAirRange(), in->getGroundRange()))*scale.first);
	}
	for (auto in : pl2->unitList)
	{
        if(in->getHealth() < EPS)
        {
            continue;
        }
		std::pair<double,double> tmp;
        tmp = std::make_pair(in->getPos().x*scale.first, in->getPos().y*scale.second);
		coord2.push_back(tmp);
		sizes.push_back((in->getSize())*scale.first);
		ranges.push_back((std::max(in->getAirRange(), in->getGroundRange()))*scale.first);
	}

	mPlayGround->setUnitsPlayer(coord1, coord2, sizes, ranges);

	//the stats have to go to the playertoolbar
	std::vector<std::string> names; //different units
	std::vector<std::string> stats; //stat-names for all units
	std::vector<std::string> values; //corresponding values
	int i = 0;
	stats.push_back("Health");
	stats.push_back("GDPS");
	stats.push_back("ADPS");
	for (auto in : pl1->unitList)
	{
        if(in->getHealth() < EPS)
        {
            continue;
        }
		std::string name = in->getName();
		names.push_back(name);
		i++;

		values.push_back(std::to_string(in->getHealth()));
		values.push_back(std::to_string(in->getGdps()));
		values.push_back(std::to_string(in->getAdps()));
	}
	mPlayer1->setUnitNames(names);
	mPlayer1->setUnitStats(stats);
	mPlayer1->setStatValues(values);

	i = 0;
	names.clear();
	values.clear();
	for (auto in : pl2->unitList)
	{
        if(in->getHealth() < EPS)
        {
            continue;
        }
		std::string name = in->getName();
		names.push_back(name);

		values.push_back(std::to_string(in->getHealth()));
		values.push_back(std::to_string(in->getGdps()));
		values.push_back(std::to_string(in->getAdps()));
	}
	mPlayer2->setUnitNames(names);
	mPlayer2->setUnitStats(stats);
	mPlayer2->setStatValues(values);
}

void GuiWindow::setStartFunc(bool (*func)(size_t const b))
{
	startFunc = func;
}

template<>
void GuiWindow::setUnits(PlayerState<Protoss> a, PlayerState<Protoss> b)
{
	setUnits(&a,&b);
}

template<>
void GuiWindow::setUnits(PlayerState<Zerg> a, PlayerState<Zerg> b)
{
	setUnits(&a,&b);
}

template<>
void GuiWindow::setUnits(PlayerState<Terran> a, PlayerState<Terran> b)
{
	setUnits(&a,&b);
}

template<>
void GuiWindow::setUnits(PlayerState<Protoss> a, PlayerState<Zerg> b)
{
	setUnits(&a,&b);
}

template<>
void GuiWindow::setUnits(PlayerState<Protoss> a, PlayerState<Terran> b)
{
	setUnits(&a,&b);
}

template<>
void GuiWindow::setUnits(PlayerState<Terran> a, PlayerState<Protoss> b)
{
	setUnits(&a,&b);
}

template<>
void GuiWindow::setUnits(PlayerState<Terran> a, PlayerState<Zerg> b)
{
	setUnits(&a,&b);
}

template<>
void GuiWindow::setUnits(PlayerState<Zerg> a, PlayerState<Protoss> b)
{
	setUnits(&a,&b);
}

template<>
void GuiWindow::setUnits(PlayerState<Zerg> a, PlayerState<Terran> b)
{
	setUnits(&a,&b);
}

void GuiWindow::playButtonPressed()
{
	std::cout << "Play Button Pressed" << std::endl;
	if (playable==true)
        updateIntervals = 100;
}

void GuiWindow::pauseButtonPressed()
{
	std::cout << "Pause Button Pressed" << std::endl;
	if (playable==true)
		updateIntervals = 0;
}

void GuiWindow::forwardButtonPressed()
{
	std::cout << "Forward Button Pressed" << std::endl;
	if (playable==true)
        updateIntervals = 2000;
}

void GuiWindow::coordinatesPlayer1(int coord)
{
	int x = coord >> 16;
	int y = coord & 0x0000FFFF;
	mPlayer1BaseCoords = std::make_pair(x,y);
	if (mPlayer2BaseCoords.first != -1)
		playable = true;
}

void GuiWindow::coordinatesPlayer2(int coord)
{
	int x = coord >> 16;
	int y = coord & 0x0000FFFF;
	mPlayer2BaseCoords = std::make_pair(x,y);
	if (mPlayer1BaseCoords.first != -1)
		playable = true;
}

void GuiWindow::showUnits()
{
	mPlayer1->showUnits();
	mPlayer2->showUnits();
}

void GuiWindow::initGuiElements()
{
	std::pair<double,double> size = mPlayGround->getPlayGroundSize();
	mPlayGround->setFixedSize(size.first+10,size.second+10);
	mPlayer1->setHeader("Player1");
	mPlayer2->setHeader("Player2");

	mPlayButton = new QPushButton("Play", this);
	mPauseButton = new QPushButton("Pause", this);
	mForwardButton = new QPushButton("Forward", this);
	mResetButton = new QPushButton("Reset", this);

	mToolbarLayout->addWidget(mPlayButton);
	mToolbarLayout->addWidget(mPauseButton);
	mToolbarLayout->addWidget(mForwardButton);
	mToolbarLayout->addWidget(mResetButton);

	mCentralLayout->addWidget(mPlayGround);
	mCentralLayout->addItem(mToolbarLayout);

	mHLayout->addWidget(mPlayer1);
	mHLayout->addItem(mCentralLayout);
	mHLayout->addWidget(mPlayer2);

}

void GuiWindow::resetButtonPressed()
{
	std::pair<int,int> initPos;
	initPos = std::make_pair(-1,-1);
	playable = true;
	mPlayer1BaseCoords = initPos;
	mPlayer2BaseCoords = initPos;
}

void GuiWindow::makeConnections()
{
	//connect ResetButton with the other classes, to reset their state
	connect(mResetButton, SIGNAL(clicked()), mPlayGround, SLOT(resetButtonPressed()));
	connect(mResetButton, SIGNAL(clicked()), mPlayer1, SLOT(resetButtonPressed()));
	connect(mResetButton, SIGNAL(clicked()), mPlayer2, SLOT(resetButtonPressed()));
	connect(mResetButton, SIGNAL(clicked()), this, SLOT(resetButtonPressed()));

	//connect Play, Pause and ForwardButton
	connect(mPlayButton, SIGNAL(clicked()), this, SLOT(playButtonPressed()));
	connect(mPauseButton, SIGNAL(clicked()), this, SLOT(pauseButtonPressed()));
	connect(mForwardButton, SIGNAL(clicked()), this, SLOT(forwardButtonPressed()));

	//connect PlaygroundSignal, when coordinates for the bases are set
	connect(mPlayGround, SIGNAL(coordinatesPlayer1(int)), this, SLOT(coordinatesPlayer1(int)));
	connect(mPlayGround, SIGNAL(coordinatesPlayer2(int)), this, SLOT(coordinatesPlayer2(int)));

	//connect PlayerToolbar, when Name was doubleclicked
	connect(mPlayer1, SIGNAL(itemDoubleClicked(QTreeWidgetItem*, int)), mPlayer1, SLOT(unitMarked(QTreeWidgetItem*)));
	connect(mPlayer2, SIGNAL(itemDoubleClicked(QTreeWidgetItem*, int)), mPlayer2, SLOT(unitMarked(QTreeWidgetItem*)));

	//connect PlayGround with PlayerToolbar, when a Unit is marked
	connect(mPlayer1, SIGNAL(markedUnit(int)), mPlayGround, SLOT(colorUnit1(int)));
	connect(mPlayer2, SIGNAL(markedUnit(int)), mPlayGround, SLOT(colorUnit2(int)));

	//connect mainThread(GUI) with side threads
	connect(this, SIGNAL(crossThreadShowUnits()), this, SLOT(showUnits()));
}
