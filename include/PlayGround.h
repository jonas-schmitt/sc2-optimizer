#ifndef PLAYGROUND_H_
#define PLAYGROUND_H_

#include <QtGui>
#include <QGraphicsView>

#include <Utilities.h>
#include <iostream>

class PlayGround : public QGraphicsView
{
	Q_OBJECT

	public:
		PlayGround(QWidget* parent = 0);

		void setUnitsPlayer(std::vector<std::pair<double,double>> unitCoordinates1, std::vector<std::pair<double,double>> unitCoordinates2, std::vector<double> sizes, std::vector<double> attackRange);
		std::pair<double, double> getPlayGroundSize();

	public slots:
		void resetButtonPressed();
		void colorUnit1(int idx);
		void colorUnit2(int idx);

	protected:
		void mousePressEvent(QMouseEvent*);

	private:
		QGraphicsScene* mScene;
		QPixmap mPixmap;

		bool basePlayer1;
		bool basePlayer2;
		double sizeX;
		double sizeY;
		size_t prev1;
		size_t prev2;

		QPoint base1;
		QPoint base2;

		std::vector<std::pair<double,double>> unitsPlayer1;
		std::vector<std::pair<double,double>> unitsPlayer2;

		//must be filled by playerToolbar
		std::vector<bool> idx1;
		std::vector<bool> idx2;

	signals:
		void coordinatesPlayer1(int);
		void coordinatesPlayer2(int);
};

#endif
