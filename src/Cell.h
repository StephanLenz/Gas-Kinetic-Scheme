
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
	// Cell center
	double centerX;
	double centerY;

	// Cell size
	double dx;
	double dy;

    FluidParameter fluidParam;

	// links to interfaces
	//    -----------
	//    |    3    |
	//    | 0     2 |
	//    |    1    |
	//    -----------
	Interface** InterfaceList;

	// Primary Variables
	double prim[4];

	// Conseved Variables
	double cons[4];

    ConservedVariable residual;

    // Boundary Cell
    BoundaryCondition* BoundaryContitionPointer;

    InterfaceType interfaceType;

public:
	Cell();
	Cell(InterfaceType interfacetype, double centerX, double centerY, double dx, double dy, BoundaryCondition* BC, FluidParameter fluidParam);

	~Cell();

	void addInterface(Interface* newInterface, int direction);

	void update(double dt);

    void applyBoundaryCondition();

	void setValues(double rho, double u, double v, double T);

    void computePrim();

    void computeCons();

    double getLocalTimestep();

	float2 getCenter();

    PrimaryVariable getPrim();

    ConservedVariable getCons();

    ConservedVariable getLocalResidual();

    float2 getDx();

    Cell* getNeighborCell(int i);

    Cell* getOpposingCell(Interface* askingInterface);

    Cell* findNeighborInDomain();

    bool isGhostCell();

	string toString();

    string valuesToString();

	string writeNodes();
};

#endif