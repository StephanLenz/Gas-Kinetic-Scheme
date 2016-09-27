

#include "Interface.h"
#include "CompressibleInterface.h"
#include <sstream>

CompressibleInterface::CompressibleInterface()
{
}

CompressibleInterface::CompressibleInterface(Cell* negCell, Cell* posCell, bool negAdd, bool posAdd, float2** nodes, FluidParameter fluidParam, BoundaryCondition* BC)
    : Interface(negCell, posCell, negAdd, posAdd, nodes, fluidParam, BC)
{
}

CompressibleInterface::~CompressibleInterface()
{
}

void CompressibleInterface::computeTimeDerivative(double * prim, double * MomentU, double * MomentV, double * MomentXi,
                                      double* a, double* b, double * timeGrad)
{
    // ========================================================================
    timeGrad[0] = a[0] * MomentU[1] * MomentV[0]
                + a[1] * MomentU[2] * MomentV[0]
                + a[2] * MomentU[1] * MomentV[1]
                + a[3] * 0.5 * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                + b[0] * MomentU[0] * MomentV[1]
                + b[1] * MomentU[1] * MomentV[1]
                + b[2] * MomentU[0] * MomentV[2]
                + b[3] * 0.5 * ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ 2.0 * prim[3] * ( MomentU[0]*MomentV[0] * prim[1] - MomentU[1]*MomentV[0] ) * this->fluidParam.Force.x
                //+ 2.0 * prim[3] * ( MomentU[0]*MomentV[0] * prim[2] - MomentU[0]*MomentV[1] ) * this->fluidParam.Force.y
                ;
    // ========================================================================
    
    // ========================================================================
    timeGrad[1] = a[0] * MomentU[2] * MomentV[0]
                + a[1] * MomentU[3] * MomentV[0]
                + a[2] * MomentU[2] * MomentV[1]
                + a[3] * 0.5 * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                + b[0] * MomentU[1] * MomentV[1]
                + b[1] * MomentU[2] * MomentV[1]
                + b[2] * MomentU[1] * MomentV[2]
                + b[3] * 0.5 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ 2.0 * prim[3] * ( MomentU[1]*MomentV[0] * prim[1] - MomentU[2]*MomentV[0] ) * this->fluidParam.Force.x
                //+ 2.0 * prim[3] * ( MomentU[1]*MomentV[0] * prim[2] - MomentU[1]*MomentV[1] ) * this->fluidParam.Force.y
                ;
    // ========================================================================
    
    // ========================================================================
    timeGrad[2] = a[0] * MomentU[1] * MomentV[1]
                + a[1] * MomentU[2] * MomentV[1]
                + a[2] * MomentU[1] * MomentV[2]
                + a[3] * 0.5 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + b[0] * MomentU[0] * MomentV[2]
                + b[1] * MomentU[1] * MomentV[2]
                + b[2] * MomentU[0] * MomentV[3]
                + b[3] * 0.5 * ( MomentU[2] * MomentV[2] + MomentU[0] * MomentV[4] + MomentU[0] * MomentV[2] * MomentXi[2] )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ 2.0 * prim[3] * ( MomentU[0]*MomentV[1] * prim[1] - MomentU[1]*MomentV[1] ) * this->fluidParam.Force.x
                //+ 2.0 * prim[3] * ( MomentU[0]*MomentV[1] * prim[2] - MomentU[0]*MomentV[2] ) * this->fluidParam.Force.y
                ;
    // ========================================================================
    
    // ========================================================================
    timeGrad[3] = a[0] * 0.50 * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                + a[1] * 0.50 * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                + a[2] * 0.50 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + a[3] * 0.25 * ( MomentU[5] + MomentU[1]* ( MomentV[4] + MomentXi[4] )
                                + 2.0 * MomentU[3] * MomentV[2]
                                + 2.0 * MomentU[3] * MomentXi[2]
                                + 2.0 * MomentU[1] * MomentV[2] * MomentXi[2]
                                )
                + b[0] * 0.50 * ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] )
                + b[1] * 0.50 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + b[2] * 0.50 * ( MomentU[2] * MomentV[2] + MomentU[0] * MomentV[4] + MomentU[0] * MomentV[2] * MomentXi[2] )
                + b[3] * 0.25 * ( MomentV[5] + MomentV[1] * ( MomentU[4] + MomentXi[4] )
                                + 2.0 * MomentU[2] * MomentV[3]
                                + 2.0 * MomentU[2] * MomentV[1] * MomentXi[2]
                                + 2.0 * MomentV[3] * MomentXi[2]
                                )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ prim[3] * ( ( MomentU[2] + MomentV[2] + MomentXi[2] ) * prim[1] 
                //            - ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                //            ) * this->fluidParam.Force.x
                //+ prim[3] * ( ( MomentU[2] + MomentV[2] + MomentXi[2] ) * prim[2] 
                //            - ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] ) 
                //            ) * this->fluidParam.Force.y
                ;
    // ========================================================================

    timeGrad[0] *= -1.0;
    timeGrad[1] *= -1.0;
    timeGrad[2] *= -1.0;
    timeGrad[3] *= -1.0;

}

void CompressibleInterface::assembleFlux(double * MomentU, double * MomentV, double * MomentXi, double * a, double * b, double * A, double * timeCoefficients, double* prim, double tau)
{
    double Flux_1[4];
    double Flux_2[4];
    double Flux_3[4];

    int u = 1;
    int v = 0;
    
    // ========================================================================
    Flux_1[0] = MomentU[0+u] * MomentV[0+v];
    Flux_1[1] = MomentU[1+u] * MomentV[0+v];
    Flux_1[2] = MomentU[0+u] * MomentV[1+v];
    Flux_1[3] = 0.5 * ( MomentU[2+u] * MomentV[0+v]
                      + MomentU[0+u] * MomentV[2+v] 
                      + MomentU[0+u] * MomentV[0+v] * MomentXi[2] );
    // ========================================================================
    
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================

    // ========================================================================
    Flux_2[0] = ( a[0] * MomentU[1+u] * MomentV[0+v]
                + a[1] * MomentU[2+u] * MomentV[0+v]
                + a[2] * MomentU[1+u] * MomentV[1+v]
                + a[3] * 0.5 * ( MomentU[3+u] * MomentV[0+v]
                               + MomentU[1+u] * MomentV[2+v] 
                               + MomentU[1+u] * MomentV[0+v] * MomentXi[2] )
                + b[0] * MomentU[0+u] * MomentV[1+v]
                + b[1] * MomentU[1+u] * MomentV[1+v]
                + b[2] * MomentU[0+u] * MomentV[2+v]
                + b[3] * 0.5 * ( MomentU[2+u] * MomentV[1+v] 
                               + MomentU[0+u] * MomentV[3+v] 
                               + MomentU[0+u] * MomentV[1+v] * MomentXi[2] )
                )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ 2.0 * prim[3] * ( MomentU[0+u]*MomentV[0+v] * prim[1] - MomentU[1+u]*MomentV[0+v] ) * this->fluidParam.Force.x
                //+ 2.0 * prim[3] * ( MomentU[0+u]*MomentV[0+v] * prim[2] - MomentU[0+u]*MomentV[1+v] ) * this->fluidParam.Force.y
                ;
    // ========================================================================
    
    // ========================================================================
    Flux_2[1] = ( a[0] * MomentU[2+u] * MomentV[0+v]
                + a[1] * MomentU[3+u] * MomentV[0+v]
                + a[2] * MomentU[2+u] * MomentV[1+v]
                + a[3] * 0.5 * ( MomentU[4+u] * MomentV[0+v] 
                               + MomentU[2+u] * MomentV[2+v] 
                               + MomentU[2+u] * MomentV[0+v] * MomentXi[2] )
                + b[0] * MomentU[1+u] * MomentV[1+v]
                + b[1] * MomentU[2+u] * MomentV[1+v]
                + b[2] * MomentU[1+u] * MomentV[2+v]
                + b[3] * 0.5 * ( MomentU[3+u] * MomentV[1+v] 
                               + MomentU[1+u] * MomentV[3+v] 
                               + MomentU[1+u] * MomentV[1+v] * MomentXi[2] )
                )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ 2.0 * prim[3] * ( MomentU[1+u]*MomentV[0+v] * prim[1] - MomentU[2+u]*MomentV[0+v] ) * this->fluidParam.Force.x
                //+ 2.0 * prim[3] * ( MomentU[1+u]*MomentV[0+v] * prim[2] - MomentU[1+u]*MomentV[1+v] ) * this->fluidParam.Force.y
                ;
    // ========================================================================
    
    // ========================================================================
    Flux_2[2] = ( a[0] * MomentU[1+u] * MomentV[1+v]
                + a[1] * MomentU[2+u] * MomentV[1+v]
                + a[2] * MomentU[1+u] * MomentV[2+v]
                + a[3] * 0.5 * ( MomentU[3+u] * MomentV[1+v]
                               + MomentU[1+u] * MomentV[3+v]
                               + MomentU[1+u] * MomentV[1+v] * MomentXi[2] )
                + b[0] * MomentU[0+u] * MomentV[2+v]
                + b[1] * MomentU[1+u] * MomentV[2+v]
                + b[2] * MomentU[0+u] * MomentV[3+v]
                + b[3] * 0.5 * ( MomentU[2+u] * MomentV[2+v]
                               + MomentU[0+u] * MomentV[4+v]
                               + MomentU[0+u] * MomentV[2+v] * MomentXi[2] )
                )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ 2.0 * prim[3] * ( MomentU[0+u]*MomentV[1+v] * prim[1] - MomentU[1+u]*MomentV[1+v] ) * this->fluidParam.Force.x
                //+ 2.0 * prim[3] * ( MomentU[0+u]*MomentV[1+v] * prim[2] - MomentU[0+u]*MomentV[2+v] ) * this->fluidParam.Force.y
                ;
    // ========================================================================
    
    // ========================================================================
    Flux_2[3] = 0.5 * ( a[0] * ( MomentU[3+u] * MomentV[0+v] 
                               + MomentU[1+u] * MomentV[2+v]
                               + MomentU[1+u] * MomentV[0+v] * MomentXi[2] )
                      + a[1] * ( MomentU[4+u] * MomentV[0+v]
                               + MomentU[2+u] * MomentV[2+v]
                               + MomentU[2+u] * MomentV[0+v] * MomentXi[2] )
                      + a[2] * ( MomentU[3+u] * MomentV[1+v]
                               + MomentU[1+u] * MomentV[3+v]
                               + MomentU[1+u] * MomentV[1+v] * MomentXi[2] )
                      + a[3] * ( 0.5 * ( MomentU[5+u] * MomentV[0+v]
                                       + MomentU[1+u] * MomentV[4+v]
                                       + MomentU[1+u] * MomentV[0+v] * MomentXi[4] )
                               +       ( MomentU[3+u] * MomentV[2+v]
                                       + MomentU[3+u] * MomentV[0+v] * MomentXi[2]
                                       + MomentU[1+u] * MomentV[2+v] * MomentXi[2] ) 
                               )
                      + b[0] * ( MomentU[2+u] * MomentV[1+v] 
                               + MomentU[0+u] * MomentV[3+v]
                               + MomentU[0+u] * MomentV[1+v] * MomentXi[2] )
                      + b[1] * ( MomentU[3+u] * MomentV[1+v]
                               + MomentU[1+u] * MomentV[3+v]
                               + MomentU[1+u] * MomentV[1+v] * MomentXi[2] )
                      + b[2] * ( MomentU[2+u] * MomentV[2+v]
                               + MomentU[0+u] * MomentV[4+v]
                               + MomentU[0+u] * MomentV[2+v] * MomentXi[2] )
                      + b[3] * ( 0.5 * ( MomentU[4+u] * MomentV[1+v] 
                                       + MomentU[0+u] * MomentV[5+v]
                                       + MomentU[0+u] * MomentV[1+v] * MomentXi[4] )
                               +       ( MomentU[2+u] * MomentV[3+v]
                                       + MomentU[2+u] * MomentV[1+v] * MomentXi[2]
                                       + MomentU[0+u] * MomentV[3+v] * MomentXi[2] )
                               )
                      )
    // this part comes from the inclusion of the forcing into the flux computation
                    //+ prim[3] * ( ( MomentU[2+u] * MomentV[0+v] + MomentU[0+u] * MomentV[2+v] + MomentU[0+u] * MomentV[0+v] * MomentXi[2] ) * prim[1] 
                    //            - ( MomentU[3+u] * MomentV[0+v] + MomentU[1+u] * MomentV[2+v] + MomentU[1+u] * MomentV[0+v] * MomentXi[2] )
                    //            ) * this->fluidParam.Force.x
                    //+ prim[3] * ( ( MomentU[2+u] * MomentV[0+v] + MomentU[0+u] * MomentV[2+v] + MomentU[0+u] * MomentV[0+v] * MomentXi[2] ) * prim[2] 
                    //            - ( MomentU[2+u] * MomentV[1+v] + MomentU[0+u] * MomentV[3+v] + MomentU[0+u] * MomentV[1+v] * MomentXi[2] ) 
                    //            ) * this->fluidParam.Force.y
                    ;
    // ========================================================================
    
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================

    // ========================================================================
    Flux_3[0] = ( A[0] * MomentU[0+u] * MomentV[0+v]
                + A[1] * MomentU[1+u] * MomentV[0+v]
                + A[2] * MomentU[0+u] * MomentV[1+v]
                + A[3] * 0.5 * ( MomentU[2+u]*MomentV[0+v]
                               + MomentU[0+u]*MomentV[2+v]
                               + MomentU[0+u]*MomentV[0+v]*MomentXi[2] )
                );
    // ========================================================================
    
    // ========================================================================
    Flux_3[1] = ( A[0] * MomentU[1+u] * MomentV[0+v]
                + A[1] * MomentU[2+u] * MomentV[0+v]
                + A[2] * MomentU[1+u] * MomentV[1+v]
                + A[3] * 0.5 * ( MomentU[3+u]*MomentV[0+v]
                               + MomentU[1+u]*MomentV[2+v]
                               + MomentU[1+u]*MomentV[0+v]*MomentXi[2] )
                );
    // ========================================================================
    
    // ========================================================================
    Flux_3[2] = ( A[0] * MomentU[0+u] * MomentV[1+v]
                + A[1] * MomentU[1+u] * MomentV[1+v]
                + A[2] * MomentU[0+u] * MomentV[2+v]
                + A[3] * 0.5 * ( MomentU[2+u]*MomentV[1+v]
                               + MomentU[0+u]*MomentV[3+v]
                               + MomentU[0+u]*MomentV[1+v]*MomentXi[2] )
                );
    // ========================================================================
    
    // ========================================================================
    Flux_3[3] = 0.5 * ( A[0] * ( MomentU[2+u] * MomentV[0+v]
                               + MomentU[0+u] * MomentV[2+v]
                               + MomentU[0+u] * MomentV[0+v] * MomentXi[2] )
                      + A[1] * ( MomentU[3+u] * MomentV[0+v]
                               + MomentU[1+u] * MomentV[2+v]
                               + MomentU[1+u] * MomentV[0+v] * MomentXi[2] )
                      + A[2] * ( MomentU[2+u] * MomentV[1+v]
                               + MomentU[0+u] * MomentV[3+v]
                               + MomentU[0+u] * MomentV[1+v] * MomentXi[2] )
                      + A[3] * ( 0.5 * ( MomentU[4+u] * MomentV[0+v]
                                       + MomentU[0+u] * MomentV[4+v]
                                       + MomentU[0+u] * MomentV[0+v] * MomentXi[4] )
                               +       ( MomentU[2+u] * MomentV[2+v]
                                       + MomentU[2+u] * MomentV[0+v] * MomentXi[2]
                                       + MomentU[0+u] * MomentV[2+v] * MomentXi[2] ) 
                               )
                      );
    // ========================================================================
    
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================

    // ========================================================================
    for ( int i = 0; i < 4; i++ )
    {
        this->timeIntegratedFlux[i]   = ( timeCoefficients[0] * Flux_1[i] + timeCoefficients[1] * Flux_2[i] + timeCoefficients[2] * Flux_3[i] ) * this->area * prim[0];

        // several test quantities
        this->timeIntegratedFlux_1[i] =   timeCoefficients[0] * Flux_1[i]                                                                       * this->area * prim[0];
        this->timeIntegratedFlux_2[i] =                                     timeCoefficients[1] * Flux_2[i]                                     * this->area * prim[0];
        this->timeIntegratedFlux_3[i] =                                                                       timeCoefficients[2] * Flux_3[i]   * this->area * prim[0];
        // The Flux density in the Flux per unit area of the interface at one instant in time
        this->FluxDensity[i] = ( Flux_1[i] - tau*( Flux_2[i] + Flux_3[i] ) ) * prim[0];
    }
    // ========================================================================

    int i = 1;
}

void CompressibleInterface::computeMicroSlope(double * prim, double * macroSlope, double * microSlope)
{
    // this method computes the micro slopes from the slopes of the conservative variables
    // the resulting microslopes contain the density, since they are computed from the slopes
    // of the conservative variables, which are rho, rhoU, rhoV and rhoE

    double A, B, C, E;

    // ========================================================================
    // this is 2 times the total energy density 2 E = 2 rhoE / rho
    E = prim[1] * prim[1] + prim[2] * prim[2] + ( this->fluidParam.K + 2.0 ) / ( 2.0*prim[3] );
    // ========================================================================

    // ========================================================================
    // the product rule of derivations is used here!
    A = 2.0*macroSlope[3] - E       * macroSlope[0];    // = 2 rho dE/dx
    B =     macroSlope[1] - prim[1] * macroSlope[0];    // =   rho dU/dx
    C =     macroSlope[2] - prim[2] * macroSlope[0];    // =   rho dV/dx
    // ========================================================================

    // compute micro slopes of primitive variables from macro slopes of conservative variables
    microSlope[3] = ( 4.0 * prim[3] * prim[3] ) / ( this->fluidParam.K + 2.0 )
                  * ( A - 2.0*prim[1] * B - 2.0*prim[2] * C );

    microSlope[2] = 2.0 * prim[3] * C - prim[2] * microSlope[3];

    microSlope[1] = 2.0 * prim[3] * B - prim[1] * microSlope[3];

    microSlope[0] = macroSlope[0] - prim[1] * microSlope[1] - prim[2] * microSlope[2] - 0.5 * E* microSlope[3];
}