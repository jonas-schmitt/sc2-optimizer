
#include "../src/InitPlayerUnits.cpp"
#include "../src/MicroSimulation.cpp"


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

