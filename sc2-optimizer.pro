######################################################################
# Automatically generated by qmake (3.0) Mi. Apr. 30 09:00:44 2014
######################################################################

TEMPLATE = app
TARGET = bin/opt
INCLUDEPATH += . include src
OBJECTS_DIR = obj
MOC_DIR = obj
QT += widgets core gui
QMAKE_CXXFLAGS += -D_GLIBCXX_USE_NANOSLEEP -std=c++11 -Wall -Wextra -O3

# Input
HEADERS += include/DataReader.h \
           include/Debug.h \
           include/GuiWindow.h \
           include/InitPlayerUnits.h \
           include/MicroSimulation.h \
           include/PlayerState.h \
	   include/PlayerToolbar.h \
	   include/PlayGround.h \
	   include/BaseSelector.h \
	   include/GuiInterface.h \
           include/Race.h \
           include/raceSelector.h \
           include/TemplateInit.h \
           include/Unit.h \
           include/UnitFactory.h \ 
	   include/UnitOptimizer.h \
	   include/UnitOptimizerBase.h \
           include/UnitGenes.h \
           include/Utilities.h
#           src/InitPlayerUnits.cpp \
#           src/MicroSimulation.cpp
SOURCES += src/DataReader.cpp \
           src/GuiWindow.cpp \
	   src/raceSelector.cpp \
           src/InitPlayerUnits.cpp \
           src/Unit.cpp \
           src/MicroSimulation.cpp \
	   src/GuiInterface.cpp \
	   src/PlayerToolbar.cpp \
           src/PlayGround.cpp \
	   src/UnitOptimizer.cpp \
	   src/main.cpp
