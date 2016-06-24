

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

void IncompressibleInterface::computeFlux(double dt)
{
    if ( this->isBoundaryInterface() )
    {
        this->computeBoundaryFlux(dt);
    }
    else
    {
        this->computeInternalFlux(dt);
    }

}

void IncompressibleInterface::computeInternalFlux(double dt)
{
    const int NUMBER_OF_MOMENTS = 7;

    double prim[4];
    //double cons[4];
    double normalGradCons[4];
    double tangentialGradCons[4];
    double timeGrad[4];

    double a[4];
    double b[4];
    double A[4];

    double MomentU[NUMBER_OF_MOMENTS];
    double MomentV[NUMBER_OF_MOMENTS];
    double MomentXi[NUMBER_OF_MOMENTS];

    // compute the length of the interface
    double dy = this->posCell->getDx().x * normal.y
        + this->posCell->getDx().y * normal.x;

    // ========================================================================
    // interpolated primary variables at the interface
    //this->interpolatePrim(prim);
    this->interpolatePrimThirdOrder(prim);

    // spacial gradients of the conservative varibles
    //this->differentiateConsNormal(normalGradCons, prim);
    this->differentiateConsNormalThirdOrder(normalGradCons, prim);
    this->differentiateConsTangential(tangentialGradCons, prim);
    // ========================================================================

    // ========================================================================
    // Formular as in the Rayleigh-Bernard-Paper (Xu, Lui, 1999)
    double tau = 2.0*prim[3] * this->fluidParam.nu;
    // ========================================================================

    // time integration Coefficients
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };

    // in case of horizontal interface (G interface), swap velocity directions
    if ( this->axis == 1 )
    {
        this->rotate(prim);
        this->rotate(normalGradCons);
        this->rotate(tangentialGradCons);
    }

    // ========================================================================
    // spacial micro slopes a = a1 + a2 u + a3 v
    //                      b = b1 + b2 u + b3 v
    // The microslopes contain the density (in opposition to Weidong Li's Code)
    this->computeMicroSlope(prim, normalGradCons,     a);
    this->computeMicroSlope(prim, tangentialGradCons, b);
    // ========================================================================

    // This Block can turn off the usage of tangential Derivatives
    // ========================================================================
    // ========================================================================
    // ========================================================================
    //for ( int i = 0; i < 4; i++ )
    //    b[i] = 0.0;
    // ========================================================================
    // ========================================================================
    // ========================================================================


    // ========================================================================
    // The Moments are computed without the density
    // Therefore their dimensions are powers of velocity
    this->computeMoments(prim, MomentU, MomentV, MomentXi, NUMBER_OF_MOMENTS);
    // ========================================================================

    // ========================================================================
    // temporal micro slopes A = A1 + A2 u + A3 v
    // The temporal macro slopes (timeGrad) also contain the density by explicit multiplication
    this->computeTimeDerivative(prim, MomentU, MomentV, MomentXi, a, b, timeGrad);

    this->computeMicroSlope(prim, timeGrad, A);
    // ========================================================================

    // ========================================================================
    // compute mass and momentum fluxes
    this->assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients, dy, prim, tau);
    // ========================================================================

    // in case of horizontal interface (G interface), swap velocity fluxes
    if ( this->axis == 1 )
    {
        this->rotate(this->timeIntegratedFlux);
        this->rotate(this->FluxDensity);
    }

}

void IncompressibleInterface::computeBoundaryFlux(double dt)
{
    PrimaryVariable prim = this->getCellInDomain()->getPrim();

    if ( this->axis == 1 )
    {
        this->rotate((double*)&prim);
    }
    
    double distance = this->distance( this->getCellInDomain()->getCenter() );

    // compute the length of the interface
    double dy = this->getCellInDomain()->getDx().x * normal.y
              + this->getCellInDomain()->getDx().y * normal.x;

    //ConservedVariable FluxDensity = this->BoundaryConditionPointer->computeBoundaryInterfaceFlux(prim, dx, this->fluidParam.nu);
    ConservedVariable FluxDensity;
    double sign = 1.0;

    if ( posCell == NULL )
        sign = -1.0;

    FluxDensity.rho  = 0.0;
    FluxDensity.rhoU = prim.rho / ( 2.0 * prim.L );
    FluxDensity.rhoV = - sign * this->fluidParam.nu * prim.rho * ( prim.V - this->BoundaryConditionPointer->getWallVelocity() ) / distance;
    FluxDensity.rhoE = 0.0;
    
    this->timeIntegratedFlux[0] = FluxDensity.rho  * dt * dy;
    this->timeIntegratedFlux[1] = FluxDensity.rhoU * dt * dy;
    this->timeIntegratedFlux[2] = FluxDensity.rhoV * dt * dy;
    this->timeIntegratedFlux[3] = FluxDensity.rhoE * dt * dy;

    this->FluxDensity[0] = FluxDensity.rho;
    this->FluxDensity[1] = FluxDensity.rhoU;
    this->FluxDensity[2] = FluxDensity.rhoV;
    this->FluxDensity[3] = FluxDensity.rhoE;

    if ( this->axis == 1 )
    {
        this->rotate(this->FluxDensity);
        this->rotate(this->timeIntegratedFlux);
    }

}

void IncompressibleInterface::computeTimeDerivative(double * prim, double * MomentU, double * MomentV, double * MomentXi,
                                      double* a, double* b, double * timeGrad)
{

    timeGrad[0] = a[0] * MomentU[1] * MomentV[0]
                + a[1] * MomentU[2] * MomentV[0]
                + a[2] * MomentU[1] * MomentV[1]
                + a[3] * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                + b[0] * MomentU[0] * MomentV[1]
                + b[1] * MomentU[1] * MomentV[1]
                + b[2] * MomentU[0] * MomentV[2]
                + b[3] * ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] ) ;

    timeGrad[1] = a[0] * MomentU[2] * MomentV[0]
                + a[1] * MomentU[3] * MomentV[0]
                + a[2] * MomentU[2] * MomentV[1]
                + a[3] * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                + b[0] * MomentU[1] * MomentV[1]
                + b[1] * MomentU[2] * MomentV[1]
                + b[2] * MomentU[1] * MomentV[2]
                + b[3] * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] );

    timeGrad[2] = a[0] * MomentU[1] * MomentV[1]
                + a[1] * MomentU[2] * MomentV[1]
                + a[2] * MomentU[1] * MomentV[2]
                + a[3] * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + b[0] * MomentU[0] * MomentV[2]
                + b[1] * MomentU[1] * MomentV[2]
                + b[2] * MomentU[0] * MomentV[3]
                + b[3] * ( MomentU[2] * MomentV[2] + MomentU[0] * MomentV[4] + MomentU[0] * MomentV[2] * MomentXi[2] );

    timeGrad[3] = a[0] * 0.50 * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                + a[1] * 0.50 * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                + a[2] * 0.50 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + a[3] * 0.25 * ( MomentU[5] + MomentU[1]* ( MomentV[4] + MomentXi[4] )
                                + 2.0 * MomentU[3] * MomentV[2]
                                + 2.0 * MomentU[3] * MomentXi[2]
                                + 2.0 * MomentU[1] * MomentV[2] * MomentXi[2] )
                + b[0] * 0.50 * ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] )
                + b[1] * 0.50 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + b[2] * 0.50 * ( MomentU[2] * MomentV[2] + MomentU[0] * MomentV[4] + MomentU[0] * MomentV[2] * MomentXi[2] )
                + b[3] * 0.25 * ( MomentV[5] + MomentV[1] * ( MomentU[4] + MomentXi[4] )
                                + 2.0 * MomentU[2] * MomentV[3]
                                + 2.0 * MomentU[2] * MomentV[1] * MomentXi[2]
                                + 2.0 * MomentV[3] * MomentXi[2] );

    timeGrad[0] *= -1.0;
    timeGrad[1] *= -1.0;
    timeGrad[2] *= -1.0;
    timeGrad[3] *= -1.0;

    // The above computed Moments do not contain the density.
    // Therefore the density is applied seperately.

    //timeGrad[0] *= -prim[0];
    //timeGrad[1] *= -prim[0];
    //timeGrad[2] *= -prim[0];
    //timeGrad[3] *= -prim[0];
}

void IncompressibleInterface::assembleFlux(double * MomentU, double * MomentV, double * MomentXi, double * a, double * b, double * A, double * timeCoefficients, double dy, double* prim, double tau)
{
    double Flux_1[4];   // this part does not contain the density (dimension of velocity powers)
    double Flux_2[4];   // this part contains the density by the micro slopes a and b
    double Flux_3[4];   // this part contains the density by the micro slope A
    
    // ========================================================================
    Flux_1[0] = MomentU[1];
    Flux_1[1] = MomentU[2];
    Flux_1[2] = MomentU[1] * MomentV[1];
    Flux_1[3] = 0.5 * (MomentU[3] + MomentU[1] * MomentV[2] + MomentU[1] * MomentXi[2]);
    // ========================================================================

    // ========================================================================
    Flux_2[0] = ( a[0] * MomentU[2] 
                + a[1] * MomentU[3]
                + a[2] * MomentU[2] * MomentV[1]
                + a[3] * 0.5 * ( MomentU[4] + MomentU[2]*MomentV[2] + MomentU[2]*MomentXi[2] )
                + b[0] * MomentU[1] * MomentV[1]
                + b[1] * MomentU[2] * MomentV[1]
                + b[2] * MomentU[1] * MomentV[2]
                + b[3] * 0.5 * ( MomentU[3]*MomentV[1] + MomentU[1]*MomentV[3] + MomentU[1]*MomentV[1]*MomentXi[2] )
                );
    Flux_2[1] = ( a[0] * MomentU[3] 
                + a[1] * MomentU[4]
                + a[2] * MomentU[3] * MomentV[1]
                + a[3] * 0.5 * ( MomentU[5] + MomentU[3]*MomentV[2] + MomentU[3]*MomentXi[2] )
                + b[0] * MomentU[2] * MomentV[1]
                + b[1] * MomentU[3] * MomentV[1]
                + b[2] * MomentU[2] * MomentV[2]
                + b[3] * 0.5 * ( MomentU[4]*MomentV[1] + MomentU[2]*MomentV[3] + MomentU[2]*MomentV[1]*MomentXi[2] )
                );
    Flux_2[2] = ( a[0] * MomentU[2] * MomentV[1]
                + a[1] * MomentU[3] * MomentV[1]
                + a[2] * MomentU[2] * MomentV[2]
                + a[3] * 0.5 * ( MomentU[4]*MomentV[1] + MomentU[2]*MomentV[3] + MomentU[2]*MomentV[1]*MomentXi[2] )
                + b[0] * MomentU[1] * MomentV[2]
                + b[1] * MomentU[2] * MomentV[2]
                + b[2] * MomentU[1] * MomentV[3]
                + b[3] * 0.5 * ( MomentU[3]*MomentV[2] + MomentU[1]*MomentV[4] + MomentU[1]*MomentV[2]*MomentXi[2] )
                );
    Flux_2[3] = 0.5 * ( a[0] * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                      + a[1] * ( MomentU[5] * MomentV[0] + MomentU[3] * MomentV[2] + MomentU[3] * MomentV[0] * MomentXi[2] )
                      + a[2] * ( MomentU[4] * MomentV[1] + MomentU[2] * MomentV[3] + MomentU[2] * MomentV[1] * MomentXi[2] )
                      + a[3] * ( 0.5 * ( MomentU[6] * MomentV[0] + MomentU[2] * MomentV[4] + MomentU[2] * MomentV[0] * MomentXi[4] )
                               +       ( MomentU[4] * MomentV[2] + MomentU[4] * MomentV[0] * MomentXi[2] + MomentU[2] * MomentV[2] * MomentXi[2] ) )
                      + b[0] * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                      + b[1] * ( MomentU[4] * MomentV[1] + MomentU[2] * MomentV[3] + MomentU[2] * MomentV[1] * MomentXi[2] )
                      + b[2] * ( MomentU[3] * MomentV[2] + MomentU[1] * MomentV[4] + MomentU[1] * MomentV[2] * MomentXi[2] )
                      + b[3] * ( 0.5 * ( MomentU[5] * MomentV[1] + MomentU[1] * MomentV[5] + MomentU[1] * MomentV[1] * MomentXi[4] )
                               +       ( MomentU[3] * MomentV[3] + MomentU[3] * MomentV[1] * MomentXi[2] + MomentU[1] * MomentV[3] * MomentXi[2] ) )
                      );
    // ========================================================================

    // ========================================================================
    Flux_3[0] = ( A[0] * MomentU[1] * MomentV[0]
                + A[1] * MomentU[2] * MomentV[0]
                + A[2] * MomentU[1] * MomentV[1]
                + A[3] * 0.5 * ( MomentU[3]*MomentV[0] + MomentU[1]*MomentV[2] + MomentU[1]*MomentV[0]*MomentXi[2] )
                );
    Flux_3[1] = ( A[0] * MomentU[2] * MomentV[0]
                + A[1] * MomentU[3] * MomentV[0]
                + A[2] * MomentU[2] * MomentV[1]
                + A[3] * 0.5 * ( MomentU[4]*MomentV[0] + MomentU[2]*MomentV[2] + MomentU[2]*MomentV[0]*MomentXi[2] )
                );
    Flux_3[2] = ( A[0] * MomentU[1] * MomentV[1]
                + A[1] * MomentU[2] * MomentV[1]
                + A[2] * MomentU[1] * MomentV[2]
                + A[3] * 0.5 * ( MomentU[3]*MomentV[1] + MomentU[1]*MomentV[3] + MomentU[1]*MomentV[1]*MomentXi[2] )
                );
    Flux_3[3] = 0.5 * ( A[0] * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                      + A[1] * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                      + A[2] * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                      + A[3] * ( 0.5 * ( MomentU[5] * MomentV[0] + MomentU[1] * MomentV[4] + MomentU[1] * MomentV[0] * MomentXi[4] )
                               +       ( MomentU[3] * MomentV[2] + MomentU[3] * MomentXi[2] + MomentU[1] * MomentV[2] * MomentXi[2] ) )
                      );
    // ========================================================================

    // ========================================================================
    for ( int i = 0; i < 4; i++ )
    {
        // Flux_2 and Flux_3 alrdy contain the density implicitly by the micro slopes a, b and A
        // Flux_1 depends only on the moments, in which the density is cancelt out.#
        // Therefore Flux_1 must also be multiplied with the density
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

    double A, B, C, E;

    // ========================================================================
    // this is 2 times the total energy density 2 E = 2 rhoE / rho
    E = prim[1] * prim[1] + prim[2] * prim[2] + ( this->fluidParam.K + 2.0 ) / ( 2.0*prim[3] );
    // ========================================================================

    // ========================================================================
    // the product rule of derivations is used here!
    A = 2.0*macroSlope[3] - E       * macroSlope[0];    // = 2 rho dE/dx
    B = macroSlope[1] - prim[1] * macroSlope[0];    // =   rho dU/dx
    C = macroSlope[2] - prim[2] * macroSlope[0];    // =   rho dV/dx
                                                    // ========================================================================

                                                    // compute micro slopes of primitive variables from macro slopes of conservatice variables
    microSlope[3] = ( 4.0 * prim[3] * prim[3] ) / ( this->fluidParam.K + 2.0 )
        * ( A - 2.0*prim[1] * B - 2.0*prim[2] * C );

    microSlope[2] = 2.0 * prim[3] * C - prim[2] * microSlope[3];

    microSlope[1] = 2.0 * prim[3] * B - prim[1] * microSlope[3];

    microSlope[0] = macroSlope[0] - prim[1] * microSlope[1] - prim[2] * microSlope[2] - 0.5 * E* microSlope[3];
}