#include "PlayGround.h"

PlayGround::PlayGround(QWidget *parent) : QGraphicsView(parent)
{
	mScene = new QGraphicsScene(this);
	this->setScene(mScene);

	sizeX = 750;
	sizeY = 750;
	mScene->setSceneRect(QRectF(0.0,0.0,sizeX-5,sizeY-5));

	mPixmap = QPixmap(sizeX, sizeY);
	mPixmap.fill(Qt::lightGray);
	mScene->addPixmap(mPixmap);

	basePlayer1 = false;
	basePlayer2 = false;

	for (unsigned int i = 0; i < 50; ++i)
	{
		idx1.push_back(false);
		idx2.push_back(false);
	}
	prev1 = 1000;
	prev2 = 1000;
}

void PlayGround::setUnitsPlayer(std::vector<std::pair<double,double>> unitCoordinates1, std::vector<std::pair<double,double>> unitCoordinates2, std::vector<double> sizes, std::vector<double> attackRange)
{
	//reset scene and paint bases at first
	mPixmap.fill(Qt::lightGray);
	QPainter painter(&mPixmap);
	painter.setPen(Qt::red);
	painter.setBrush(Qt::red);
	painter.drawEllipse(base1, 15, 15);
	painter.setPen(Qt::blue);
	painter.setBrush(Qt::blue);
	painter.drawEllipse(base2, 15, 15);

	unitsPlayer1 = unitCoordinates1;
	unitsPlayer2 = unitCoordinates2;

	if (prev1 > unitsPlayer1.size())
	{
		prev1 = unitsPlayer1.size();
		for (unsigned int i = 0; i < 50; i++)
		{
			idx1.at(i) = false;
		}
	}

	if (prev2 > unitsPlayer2.size())
	{
		prev2 = unitsPlayer2.size();
		for (unsigned int i = 0; i < 50; i++)
		{
			idx2.at(i) = false;
		}
	}

	//paint PlayerCoordinates
	for (unsigned int i = 0; i < unitCoordinates1.size(); ++i)
	{
		if (idx1.at(i) == true)
		{
			painter.setPen(Qt::magenta);
			painter.setBrush(Qt::magenta);
		}
		else
		{
			painter.setPen(Qt::red);
			painter.setBrush(Qt::red);
		}
		QPoint p(unitCoordinates1.at(i).first, unitCoordinates1.at(i).second);
		painter.setOpacity(0.2);
		painter.drawEllipse(p, (int)attackRange.at(i), (int)attackRange.at(i));
		painter.setOpacity(1.0);
		painter.drawEllipse(p, (int)sizes.at(i), (int)sizes.at(i));
	}

	for (unsigned int i = 0; i < unitCoordinates2.size(); ++i)
	{
		if (idx2.at(i) == true)
		{
			painter.setPen(Qt::cyan);
			painter.setBrush(Qt::cyan);
		}
		else
		{
			painter.setPen(Qt::blue);
			painter.setBrush(Qt::blue);
		}
		QPoint p(unitCoordinates2.at(i).first, unitCoordinates2.at(i).second);
		painter.setOpacity(0.2);
		painter.drawEllipse(p, (int)attackRange.at(i+unitCoordinates1.size()),(int)attackRange.at(i+unitCoordinates1.size()));
		painter.setOpacity(1.0);
		painter.drawEllipse(p, (int)sizes.at(i+unitCoordinates1.size()),(int)sizes.at(i+unitCoordinates1.size()));
	}
	painter.end();
	mScene->addPixmap(mPixmap);
}

std::pair<double,double> PlayGround::getPlayGroundSize()
{
	return std::make_pair(sizeX,sizeY);
}

void PlayGround::resetButtonPressed()
{
//	std::cout << "PlayGround reseted" << std::endl;
	basePlayer1 = false;
	basePlayer2 = false;

	mPixmap.fill(Qt::lightGray);
	mScene->addPixmap(mPixmap);
}

void PlayGround::colorUnit1(int idx)
{
	QPainter painter(&mPixmap);
	if (idx1.at(idx) == true)
	{
		idx1.at(idx) = false;
		painter.setPen(Qt::red);
		painter.setBrush(Qt::red);
	}
	else
	{
		idx1.at(idx) = true;
		painter.setPen(Qt::magenta);
		painter.setBrush(Qt::magenta);
	}
	QPoint p(unitsPlayer1.at(idx).first, unitsPlayer1.at(idx).second);
	painter.drawEllipse(p, 5, 5);
	painter.end();
	mScene->addPixmap(mPixmap);
}

void PlayGround::colorUnit2(int idx)
{
	QPainter painter(&mPixmap);
	if (idx2.at(idx) == true)
	{
		idx2.at(idx) = false;
		painter.setPen(Qt::blue);
		painter.setBrush(Qt::blue);
	}
	else
	{
		idx2.at(idx) = true;
		painter.setPen(Qt::cyan);
		painter.setBrush(Qt::cyan);
	}
	QPoint p(unitsPlayer2.at(idx).first, unitsPlayer2.at(idx).second);
	painter.drawEllipse(p, 5, 5);
	painter.end();
	mScene->addPixmap(mPixmap);
}

void PlayGround::mousePressEvent(QMouseEvent* event)
{
	if (event->button() == Qt::LeftButton && basePlayer1 == false)
	{
		base1 = event->pos();
		QPainter painter(&mPixmap);
		painter.setPen(Qt::red);
		painter.setBrush(Qt::red);
		painter.drawEllipse(base1, 15, 15);
		painter.end();
		mScene->addPixmap(mPixmap);

		basePlayer1 = true;
		int coord = base1.x()*64*1024 + base1.y();
		emit(coordinatesPlayer1(coord));
	}
	else if (event->button() == Qt::RightButton && basePlayer2 == false)
	{
		base2 = event->pos();
		QPainter painter(&mPixmap);
		painter.setPen(Qt::blue);
		painter.setBrush(Qt::blue);
		painter.drawEllipse(base2, 15, 15);
		painter.end();
		mScene->addPixmap(mPixmap);

		basePlayer2 = true;
		int coord = base2.x()*64*1024 + base2.y();
		emit(coordinatesPlayer2(coord));
	}
}
