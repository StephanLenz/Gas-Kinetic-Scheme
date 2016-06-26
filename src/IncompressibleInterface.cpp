

#include "Interface.h"
#include "IncompressibleInterface.h"
#include <sstream>

IncompressibleInterface::IncompressibleInterface()
{
}

IncompressibleInterface::IncompressibleInterface(Cell* negCell, Cell* posCell, float2 center, float2 normal, FluidParameter fluidParam, InterfaceBC* BC)
    : Interface(negCell, posCell, center, normal, fluidParam, BC)
{
}

IncompressibleInterface::~IncompressibleInterface()
{
}

void IncompressibleInterface::computeTimeDerivative(double * prim, double * MomentU, double * MomentV, double * MomentXi,
                                      double* a, double* b, double * timeGrad)
{

    timeGrad[0] = a[0] * MomentU[1] * MomentV[0]
                + a[1] * MomentU[2] * MomentV[0]
                + a[2] * MomentU[1] * MomentV[1]
                + b[0] * MomentU[0] * MomentV[1]
                + b[1] * MomentU[1] * MomentV[1]
                + b[2] * MomentU[0] * MomentV[2] ;

    timeGrad[1] = a[0] * MomentU[2] * MomentV[0]
                + a[1] * MomentU[3] * MomentV[0]
                + a[2] * MomentU[2] * MomentV[1]
                + b[0] * MomentU[1] * MomentV[1]
                + b[1] * MomentU[2] * MomentV[1]
                + b[2] * MomentU[1] * MomentV[2] ;

    timeGrad[2] = a[0] * MomentU[1] * MomentV[1]
                + a[1] * MomentU[2] * MomentV[1]
                + a[2] * MomentU[1] * MomentV[2]
                + b[0] * MomentU[0] * MomentV[2]
                + b[1] * MomentU[1] * MomentV[2]
                + b[2] * MomentU[0] * MomentV[3] ;

    timeGrad[3] = 0.0;

    timeGrad[0] *= -1.0;
    timeGrad[1] *= -1.0;
    timeGrad[2] *= -1.0;
    timeGrad[3] *= -1.0;
}

void IncompressibleInterface::assembleFlux(double * MomentU, double * MomentV, double * MomentXi, double * a, double * b, double * A, double * timeCoefficients, double dy, double* prim, double tau)
{
    double Flux_1[4];
    double Flux_2[4];
    double Flux_3[4];
    if( this->axis == 0 )
    {
    // ========================================================================
        Flux_1[0] = MomentU[1];
        Flux_1[1] = MomentU[2];
        Flux_1[2] = MomentU[1] * MomentV[1];
        Flux_1[3] = 0.0;
    // ========================================================================

    // ========================================================================
        Flux_2[0] = ( a[0] * MomentU[2] * MomentV[0]
                    + a[1] * MomentU[3] * MomentV[0]
                    + a[2] * MomentU[2] * MomentV[1]
                    + b[0] * MomentU[1] * MomentV[1]
                    + b[1] * MomentU[2] * MomentV[1]
                    + b[2] * MomentU[1] * MomentV[2]
                    );
        Flux_2[1] = ( a[0] * MomentU[3] * MomentV[0]
                    + a[1] * MomentU[4] * MomentV[0]
                    + a[2] * MomentU[3] * MomentV[1]
                    + b[0] * MomentU[2] * MomentV[1]
                    + b[1] * MomentU[3] * MomentV[1]
                    + b[2] * MomentU[2] * MomentV[2]
                    );
        Flux_2[2] = ( a[0] * MomentU[2] * MomentV[1]
                    + a[1] * MomentU[3] * MomentV[1]
                    + a[2] * MomentU[2] * MomentV[2]
                    + b[0] * MomentU[1] * MomentV[2]
                    + b[1] * MomentU[2] * MomentV[2]
                    + b[2] * MomentU[1] * MomentV[3]
                    );
        Flux_2[3] = 0.0;
    // ========================================================================

    // ========================================================================
        Flux_3[0] = ( A[0] * MomentU[1] * MomentV[0]
                    + A[1] * MomentU[2] * MomentV[0]
                    + A[2] * MomentU[1] * MomentV[1]
                    );
        Flux_3[1] = ( A[0] * MomentU[2] * MomentV[0]
                    + A[1] * MomentU[3] * MomentV[0]
                    + A[2] * MomentU[2] * MomentV[1]
                    );
        Flux_3[2] = ( A[0] * MomentU[1] * MomentV[1]
                    + A[1] * MomentU[2] * MomentV[1]
                    + A[2] * MomentU[1] * MomentV[2]
                    );
        Flux_3[3] = 0.0;
    // ========================================================================
    }
    else
    {
    // ========================================================================
        Flux_1[0] = MomentU[0] * MomentV[1];
        Flux_1[1] = MomentU[1] * MomentV[1];
        Flux_1[2] = MomentU[0] * MomentV[2];
        Flux_1[3] = 0.0;
    // ========================================================================

    // ========================================================================
        Flux_2[0] = ( a[0] * MomentU[1] * MomentV[1]
                    + a[1] * MomentU[2] * MomentV[1]
                    + a[2] * MomentU[1] * MomentV[2]
                    + b[0] * MomentU[0] * MomentV[2]
                    + b[1] * MomentU[1] * MomentV[2]
                    + b[2] * MomentU[0] * MomentV[3]
                    );
        Flux_2[1] = ( a[0] * MomentU[2] * MomentV[1]
                    + a[1] * MomentU[3] * MomentV[1]
                    + a[2] * MomentU[2] * MomentV[2]
                    + b[0] * MomentU[1] * MomentV[2]
                    + b[1] * MomentU[2] * MomentV[2]
                    + b[2] * MomentU[1] * MomentV[3]
                    );
        Flux_2[2] = ( a[0] * MomentU[1] * MomentV[2]
                    + a[1] * MomentU[2] * MomentV[2]
                    + a[2] * MomentU[1] * MomentV[3]
                    + b[0] * MomentU[0] * MomentV[3]
                    + b[1] * MomentU[1] * MomentV[3]
                    + b[2] * MomentU[0] * MomentV[4]
                    );
        Flux_2[3] = 0.0;
    // ========================================================================

    // ========================================================================
        Flux_3[0] = ( A[0] * MomentU[0] * MomentV[1]
                    + A[1] * MomentU[1] * MomentV[1]
                    + A[2] * MomentU[0] * MomentV[2]
                    );
        Flux_3[1] = ( A[0] * MomentU[1] * MomentV[1]
                    + A[1] * MomentU[2] * MomentV[1]
                    + A[2] * MomentU[1] * MomentV[2]
                    );
        Flux_3[2] = ( A[0] * MomentU[0] * MomentV[2]
                    + A[1] * MomentU[1] * MomentV[2]
                    + A[2] * MomentU[0] * MomentV[3]
                    );
        Flux_3[3] = 0.0;
    // ========================================================================
    }
    // ========================================================================
    for ( int i = 0; i < 4; i++ )
    {
        this->timeIntegratedFlux[i] = ( timeCoefficients[0] * Flux_1[i] + timeCoefficients[1] * Flux_2[i] + timeCoefficients[2] * Flux_3[i] ) * dy * prim[0];
        // The Flux density in the Flux per unit area of the interface at one instant in time
        this->FluxDensity[i] = ( Flux_1[i] - tau*( Flux_2[i] + Flux_3[i] ) ) * prim[0];
    }
    // ========================================================================
}

void IncompressibleInterface::computeMicroSlope(double * prim, double * macroSlope, double * microSlope)
{
    // this method computes the micro slopes from the slopes of the conservative variables
    // the resulting microslopes contain the density, since they are computed from the slopes
    // of the conservative variables, which are rho, rhoU, rhoV and rhoE

    microSlope[3] = 0.0;

    microSlope[2] = 2.0 * prim[3] * ( macroSlope[2] - prim[2] * macroSlope[0] );

    microSlope[1] = 2.0 * prim[3] * ( macroSlope[1] - prim[1] * macroSlope[0] );

    microSlope[0] = macroSlope[0] - prim[1] * microSlope[1] - prim[2] * microSlope[2];
}