
#ifndef INTERFACE_H
#define INTERFACE_H

#include "Cell.h"
#include "InterfaceBC.h"
//class Cell;
#include <string>

using namespace std;

class Interface
{
protected:
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

    static Interface* createInterface(InterfaceType type, Cell* negCell, Cell* posCell, float2 center, float2 normal, FluidParameter fluidParam, InterfaceBC* BC);

	virtual void computeFlux(double dt);

    virtual void computeInternalFlux(double dt);

    virtual void computeBoundaryFlux(double dt);

    Cell* getNeigborCell(Cell* askingCell);
    Cell* getCellInDomain();

    ConservedVariable getTimeIntegratedFlux();
    ConservedVariable getFluxDensity();

    bool isGhostInterface();
    bool isBoundaryInterface();

	string toString();

    string writeCenter();

protected:

    void interpolatePrim(double* prim);
    void interpolatePrimThirdOrder(double* prim);

    void interpolateCons(double* prim);
    void interpolateConsThirdOrder(double* prim);

    void differentiateConsNormal(double* normalGradCons, double* prim);
    void differentiateConsNormalThirdOrder(double* normalGradCons, double* prim);

    void differentiateConsTangential(double* tangentialGradCons, double* prim);

    virtual void computeTimeDerivative(double* prim, double* MomentU, double* MomentV, double* MomentXi,
                                       double* a, double* b, double * timeGrad) = 0;

    virtual void assembleFlux(double* MomentU, double* MomentV, double* MomentXi, 
                              double* a, double* b, double* A, double* timeCoefficients,
                              double dy, double* prim, double tau) = 0;

    void rotate(double* vector);
    void cons2prim(double* prim, double*cons);
    double distance(float2 point);

    virtual void computeMicroSlope(double* prim, double* macroSlope, double* microSlope) = 0;
    void computeMoments(double* prim, double* MomentU, double* MomentV, double* MomentXi, int numberMoments);
};

#endif