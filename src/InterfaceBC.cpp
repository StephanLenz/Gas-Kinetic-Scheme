#include "InterfaceBC.h"
#include "Types.h"

InterfaceBC::InterfaceBC(double wallVelocity)
{
    this->wallVelocity = wallVelocity;
}


InterfaceBC::~InterfaceBC()
{
}

ConservedVariable InterfaceBC::computeBoundaryInterfaceFlux(PrimitiveVariable prim, double dx, double nu)
{

    ConservedVariable Flux;
    Flux.rho = 0.0;
    Flux.rhoU = 1.0 / ( 2.0 * prim.L );
    Flux.rhoV = -nu * prim.V / ( 0.5 *dx );
    Flux.rhoE = 0.0;

    return Flux;
}

double InterfaceBC::getWallVelocity()
{
    return this->wallVelocity;
}
