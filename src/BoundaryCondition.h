#include "Types.h"
#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

class BoundaryCondition
{
private:
    BoundaryConditionType type;
    PrimitiveVariable value;
public:
    BoundaryCondition();
    BoundaryCondition(  BoundaryConditionType type,
                        double rho, double U, double V, double T);
    ~BoundaryCondition();

    BoundaryConditionType getType();

    PrimitiveVariable getValue();
};

#endif

