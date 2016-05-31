
#ifndef INTERFACE_H
#define INTERFACE_H

#include "Cell.h"
#include "InterfaceBC.h"
//class Cell;
#include <string>

using namespace std;

class Interface
{
private:
	Cell* negCell;
	Cell* posCell;

    float2 center;
    float2 normal;
    int axis;

    InterfaceBC* BoundaryConditionPointer;

    FluidParameter fluidParam;

    double timeIntegratedFlux[4];
    double FluxDensity[4];
public:
	Interface();
	Interface(Cell* negCell, Cell* posCell, float2 center, float2 normal, FluidParameter fluidParam, InterfaceBC* BC);
	~Interface();

	void computeFlux(double dt);

    Cell* getNeigborCell(Cell* askingCell);
    ConservedVariable getTimeIntegratedFlux();
    ConservedVariable getFluxDensity();

    bool isGhostInterface();
    bool isBoundaryInterface();

	string toString();

    string writeCenter();

private:

    void interpolatePrim(double* prim);
    void differentiateCons(double* normalGradCons, double* tangentialGradCons, double* prim);

    void computeTimeDerivative(double* prim, double* MomentU, double* MomentV, double* MomentXi,
                               double* a, double* b, double * timeGrad);

    void assembleFlux(double* MomentU, double* MomentV, double* MomentXi, 
                      double* a, double* b, double* A, double* timeCoefficients, double dy, double* prim, double tau);

    void rotate(double* vector);

    void computeMicroSlope(double* prim, double* macroSlope, double* microSlope);
    void computeMoments(double* prim, double* MomentU, double* MomentV, double* MomentXi, int numberMoments);
};

#endif