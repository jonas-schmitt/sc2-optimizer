
#include "../src/InitPlayerUnits.cpp"
#include "../src/GuiInterface.cpp"
#include "../src/MicroSimulation.cpp"
#include "../src/UnitOptimizer.cpp"

template class InitPlayerUnits<Zerg>;
template class InitPlayerUnits<Terran>;
template class InitPlayerUnits<Protoss>;

template class MicroSimulation<Terran, Protoss>;
template class MicroSimulation<Terran, Zerg>;
template class MicroSimulation<Terran, Terran>;

template class MicroSimulation<Protoss, Zerg>;
template class MicroSimulation<Protoss, Terran>;
template class MicroSimulation<Protoss, Protoss>;

template class MicroSimulation<Zerg, Protoss>;
template class MicroSimulation<Zerg, Terran>;
template class MicroSimulation<Zerg, Zerg>;

template class GuiInterface<Terran, Protoss>;
template class GuiInterface<Terran, Zerg>;
template class GuiInterface<Terran, Terran>;

template class GuiInterface<Protoss, Zerg>;
template class GuiInterface<Protoss, Terran>;
template class GuiInterface<Protoss, Protoss>;

template class GuiInterface<Zerg, Protoss>;
template class GuiInterface<Zerg, Terran>;
template class GuiInterface<Zerg, Zerg>;

template class UnitOptimizer<Terran,Zerg>;
template class UnitOptimizer<Terran,Protoss>;

