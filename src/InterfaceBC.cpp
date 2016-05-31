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
    Flux.rho = 0.0;
    Flux.rhoU = 0.0;
    Flux.rhoV = 0.0;
    Flux.rhoE = 0.0;

    return Flux;
}
