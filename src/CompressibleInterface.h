// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      CompressibleInterface.h
//
//          This class implements the GKS specific flux computations.
//
//      Function:
//          Storage of Interface data
//          Implementation of Interface related computations
//              Flux computation 
//
// ============================================================================

#ifndef COMPRESSIBLEINTERFACE_H
#define COMPRESSIBLEINTERFACE_H

#include "Interface.h"
#include "Cell.h"
#include "InterfaceBC.h"
#include <string>

using namespace std;

class CompressibleInterface : public Interface
{
public:
	CompressibleInterface();
    CompressibleInterface(Cell* negCell, Cell* posCell, bool negAdd, bool posAdd, float2** nodes, FluidParameter fluidParam, BoundaryCondition* BC, double periodicLengthX, double periodicLengthY);
	~CompressibleInterface();

protected:

    void computeTimeDerivative(double* prim, double* MomentU, double* MomentV, double* MomentXi,
                               double* a, double* b, double * timeGrad);

    void assembleFlux(double* MomentU, double* MomentV, double* MomentXi, 
                      double* a, double* b, double* A, double* timeCoefficients, double* prim, double tau);

    void computeMicroSlope(double* prim, double* macroSlope, double* microSlope);

};

#endif