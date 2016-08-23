#include "BoundaryCondition.h"



BoundaryCondition::BoundaryCondition()
{
}

BoundaryCondition::BoundaryCondition( BoundaryConditionType type,
                                      double rho, double U, double V, double L)
{
    this->type = type;
    this->value.rho = rho;
    this->value.U = U;
    this->value.V = V;
    this->value.L = L;
}


BoundaryCondition::~BoundaryCondition()
{
}

BoundaryConditionType BoundaryCondition::getType()
{
    return this->type;
}

PrimitiveVariable BoundaryCondition::getValue()
{
    return value;
}
