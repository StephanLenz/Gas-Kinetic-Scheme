
#ifndef INCOMPRESSIBLEINTERFACE_H
#define INCOMPRESSIBLEINTERFACE_H

#include "Interface.h"
#include "Cell.h"
#include "InterfaceBC.h"
#include <string>

using namespace std;

class IncompressibleInterface : public Interface
{
public:
	IncompressibleInterface();
    IncompressibleInterface(Cell* negCell, Cell* posCell, float2 center, float2 normal, FluidParameter fluidParam, InterfaceBC* BC);
	~IncompressibleInterface();

	void computeFlux(double dt);

    void computeInternalFlux(double dt);

    void computeBoundaryFlux(double dt);

protected:

    void computeTimeDerivative(double* prim, double* MomentU, double* MomentV, double* MomentXi,
                               double* a, double* b, double * timeGrad);

    void assembleFlux(double* MomentU, double* MomentV, double* MomentXi, 
                      double* a, double* b, double* A, double* timeCoefficients, double dy, double* prim, double tau);

    void computeMicroSlope(double* prim, double* macroSlope, double* microSlope);

};

#endif