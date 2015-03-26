#ifndef PLAYERTOOLBAR_H_
#define PLAYERTOOLBAR_H_

#include <QtGui>
#include <QTreeWidget>
#include <QCheckBox>
#include <QString>

#include <string>
#include <iostream>

class PlayerToolbar : public QTreeWidget
{
	Q_OBJECT

	public:
		PlayerToolbar(QWidget *parent = 0);
		~PlayerToolbar();
		void setHeader(const QString name);
		void setUnitNames(std::vector<std::string>);
		void setUnitStats(std::vector<std::string>);
		void setStatValues(std::vector<std::string>);
		void showUnits();

	public slots:
		void resetButtonPressed();
		void unitMarked(QTreeWidgetItem* item);

	private:
		size_t lastSize;
		QList<QTreeWidgetItem*> mUnitBox;
		QList<bool> mUnitCheck;
		QList<QTreeWidgetItem*> mStatBox;

		std::vector<std::string> mNames;
		std::vector<std::string> mStats;
		std::vector<std::string> mValues;

	signals:
		void markedUnit(int);

};

#endif
