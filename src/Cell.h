
#ifndef RECTCELL2D_H
#define RECTCELL2D_H

//#include "Interface.h"
class Interface;
#include "Types.h"
#include "BoundaryCondition.h"
#include <string>

using namespace std;

class Cell
{
private:

    float2* nodes[4];

	float2 center;
    double volume;
    double minDx;

    FluidParameter fluidParam;

    int nInterfaces;
	Interface* InterfaceList[4];

	// Conseved Variables
	double cons[4];
    double cons_old[4];

    ConservedVariable residual;

    // Boundary Cell
    BoundaryCondition* BoundaryContitionPointer;

    InterfaceType interfaceType;

public:
	Cell();
	Cell(InterfaceType interfacetype, float2** nodes, BoundaryCondition* BC, FluidParameter fluidParam);

	~Cell();

	void addInterface(Interface* newInterface);

    void computeMinDx();

	void update(double dt);

    void applyBoundaryCondition();

    void applyForcing(double dt);

	void setValues(double rho, double u, double v, double T);

    void computeCons(PrimitiveVariable prim);

    double getLocalTimestep();

	float2 getCenter();
    float2 getNode(int i);

    PrimitiveVariable getPrim();

    ConservedVariable getCons();

    ConservedVariable getLocalResidual();

    Cell* getNeighborCell(int i);

    Cell* getOpposingCell(Interface* askingInterface);

    Cell* findNeighborInDomain();

    bool isGhostCell();

    double distance(float2 point);

	string toString();

    string valuesToString();

	string writeNodes();
};

#endif