
#ifndef COMPRESSIBLEINTERFACE_H
#define COMPRESSIBLEINTERFACE_H

#include "CompressibleInterface.h"
#include "Cell.h"
#include "InterfaceBC.h"
//class Cell;
#include <string>

using namespace std;

class CompressibleInterface : public Interface
{
public:
	CompressibleInterface();
    CompressibleInterface(Cell* negCell, Cell* posCell, float2 center, float2 normal, FluidParameter fluidParam, InterfaceBC* BC);
	~CompressibleInterface();

	void computeFlux(double dt);

    void computeInternalFlux(double dt);

    void computeBoundaryFlux(double dt);

protected:

    void computeTimeDerivative(double* prim, double* MomentU, double* MomentV, double* MomentXi,
                               double* a, double* b, double * timeGrad);

    void assembleFlux(double* MomentU, double* MomentV, double* MomentXi, 
                      double* a, double* b, double* A, double* timeCoefficients, double dy, double* prim, double tau);

};

#endif