#ifndef GUIWINDOW_H_
#define GUIWINDOW_H_

#include <QMainWindow>
#include <QtGui>
#include <QVBoxLayout>
#include <QPushButton>
#include <QWidget>
#include <QString>

#include <string>
#include <iostream>
#include <future>
#include <thread>

#include "PlayGround.h"
#include "PlayerToolbar.h"
#include "Race.h"
#include "PlayerState.h"
#include "Unit.h"

class GuiWindow : public QMainWindow
{
	Q_OBJECT

	public:
		GuiWindow(QWidget* parent = 0);
		~GuiWindow();

		//if the base of player 1 should be returned: true
		//else the base of player 2 is returned
		std::pair<int,int> getBaseCoordinates(bool player);
		template<typename T, typename U>
		void setUnits(T pl1, U pl2);

		void setStartFunc(bool (*func)(size_t const b));

	public slots:
		void playButtonPressed();
		void pauseButtonPressed();
		void forwardButtonPressed();
		void resetButtonPressed();
		void coordinatesPlayer1(int coord);
		void coordinatesPlayer2(int coord);
		void showUnits();

	signals:
		void crossThreadShowUnits();


	private:
		std::pair<double,double> scale;
		template<typename T, typename U>
		void setUnits(T *pl1, U *pl2);
		void scaleFactor(std::pair<double,double>size);
		std::future<void> runnerThread;
		std::thread globalThread;
		volatile int updateIntervals; //0 means 'pause'

		std::thread simThread;

		void threadWrapper();
		bool (*startFunc)(size_t const b);
		bool playable;

		void initGuiElements();
		void makeConnections();

		std::pair<int,int> mPlayer1BaseCoords;
		std::pair<int,int> mPlayer2BaseCoords;

		QWidget* mWidget;
		QHBoxLayout* mHLayout;
		QVBoxLayout* mCentralLayout;
		QHBoxLayout* mToolbarLayout;

		PlayGround* mPlayGround;
		PlayerToolbar* mPlayer1;
		PlayerToolbar* mPlayer2;

		QPushButton *mPlayButton;
		QPushButton *mPauseButton;
		QPushButton *mForwardButton;
		QPushButton *mResetButton;

};

#endif
