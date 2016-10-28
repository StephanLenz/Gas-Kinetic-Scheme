#include "Types.h"
class GKSSolver;
#include <string>
#include <vector>
#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

using namespace std;

class BoundaryCondition
{
private:
    BoundaryConditionType type;
    PrimitiveVariable value;

    vector<idType> Cell;
    vector<idType> NeighborCell;

public:
    BoundaryCondition();
    BoundaryCondition(  BoundaryConditionType type,
                        double rho, double U, double V, double T);
    ~BoundaryCondition();

    BoundaryConditionType getType();

    PrimitiveVariable getValue();

    void addCell(idType id);

    void addNeighborCell(idType id);

    void setGhostCells(GKSSolver& solver);

    string toString();
};

#endif

