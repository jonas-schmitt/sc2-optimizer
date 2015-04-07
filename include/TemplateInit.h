
#include "../src/InitPlayerUnits.cpp"
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

template class UnitOptimizer<Terran, Protoss>;
template class UnitOptimizer<Terran, Zerg>;
template class UnitOptimizer<Terran, Terran>;

template class UnitOptimizer<Protoss, Zerg>;
template class UnitOptimizer<Protoss, Terran>;
template class UnitOptimizer<Protoss, Protoss>;

template class UnitOptimizer<Zerg, Protoss>;
template class UnitOptimizer<Zerg, Terran>;
template class UnitOptimizer<Zerg, Zerg>;

