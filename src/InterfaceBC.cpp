#include "InterfaceBC.h"
#include "Types.h"

InterfaceBC::InterfaceBC(double wallVelocity)
{
    this->wallVelocity = wallVelocity;
}


InterfaceBC::~InterfaceBC()
{
}

ConservedVariable InterfaceBC::computeBoundaryInterfaceFlux(Cell * CellInDomain)
{
    ConservedVariable Flux;

    return Flux;
}
