

#ifndef INTERFACEBC_H
#define INTERFACEBC_H

#include "Types.h"
#include "Cell.h"

class InterfaceBC
{
private:

    double wallVelocity;

public:
    InterfaceBC(double wallVelocity);
    ~InterfaceBC();

    ConservedVariable computeBoundaryInterfaceFlux(PrimaryVariable prim, double dx, double nu);
};


#endif