

#include "Interface.h"
#include <sstream>

Interface::Interface()
{
}

Interface::Interface(Cell* negCell, Cell* posCell, float2 center, float2 normal, FluidParameter fluidParam, InterfaceBC* BC)
{
	this->negCell = negCell;
	this->posCell = posCell;

    this->center = center;
    this->normal = normal;

    this->BoundaryConditionPointer = BC;

    if      ( (fabs(this->normal.x - 1.0) < 1e-6) && (fabs(this->normal.y - 0.0) < 1.0e-6) )
        this->axis = 0;
    else if ( (fabs(this->normal.x - 0.0) < 1e-6) && (fabs(this->normal.y - 1.0) < 1.0e-6) )
        this->axis = 1;
    
    // links to interfaces
    //    ---------------------
    //    |    3    |    3    |
    //    | 0     2 | 0     2 |
    //    |    1    |    1    |
    //    ---------------------
    //     neg Cell   pos Cell
    //
    //    -----------
    //    |    3    |
    //    | 0     2 |   pos Cell
    //    |    1    |
    //    -----------
    //    |    3    |
    //    | 0     2 |   neg Cell
    //    |    1    |
    //    -----------

    if(this->negCell != NULL )
	    this->negCell->addInterface(this,axis+2);
    if(this->posCell != NULL )
	    this->posCell->addInterface(this,axis+0);

    this->fluidParam = fluidParam;

    for ( int i = 0; i < 4; i++ )
    {
        this->timeIntegratedFlux[i] = 0.0;
        this->FluxDensity[i] = 0.0;
    }
}

Interface::~Interface()
{
}

void Interface::computeFlux(double dt)
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

void Interface::computeInternalFlux(double dt)
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
    

    /*
    // ========================================================================
    // interpolated primary variables at the interface
    //this->interpolatePrim(cons);
    this->interpolateConsThirdOrder(cons);
    this->cons2prim(prim, cons);

    // spacial gradients of the conservative varibles
    //this->differentiateConsNormal(normalGradCons, prim);
    this->differentiateConsNormalThirdOrder(normalGradCons, prim);
    this->differentiateConsTangential(tangentialGradCons, prim);
    // ========================================================================
    */

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

void Interface::computeBoundaryFlux(double dt)
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
    FluxDensity.rhoV = - sign * this->fluidParam.nu * prim.rho * prim.V / distance;
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

Cell * Interface::getNeigborCell(Cell * askingCell)
{
    // ========================================================================
    // ============= This is only Experimental and not really correct =========
    // ========================================================================
    // In the case of an Boundary Interface, the  Cell in the Domain is 
    // virtually copied to a ghost Cell. This is wrong because this does not
    // satisfy the correct velocity at the wall. Better would be to create a
    // temporal GhostCell with the correct values!
    // ========================================================================l.
    if ( this->isBoundaryInterface() )
        return askingCell;
    // ========================================================================
    // ========================================================================
    // ========================================================================

    if (posCell == askingCell)
        return negCell;
    else
        return posCell;
}

Cell * Interface::getCellInDomain()
{
    if( !this->isBoundaryInterface() )
        return NULL;
    if ( posCell != NULL )
        return posCell;
    if ( negCell != NULL )
        return negCell;
}

ConservedVariable Interface::getTimeIntegratedFlux()
{
    ConservedVariable tmp;
    tmp.rho  = this->timeIntegratedFlux[0];
    tmp.rhoU = this->timeIntegratedFlux[1];
    tmp.rhoV = this->timeIntegratedFlux[2];
    tmp.rhoE = this->timeIntegratedFlux[3];
    return tmp;
}

ConservedVariable Interface::getFluxDensity()
{
    ConservedVariable tmp;
    tmp.rho =  this->FluxDensity[0];
    tmp.rhoU = this->FluxDensity[1];
    tmp.rhoV = this->FluxDensity[2];
    tmp.rhoE = this->FluxDensity[3];
    return tmp;
}

bool Interface::isGhostInterface()
{
    if (this->isBoundaryInterface()) return false;

    return this->posCell->isGhostCell() && this->negCell->isGhostCell();
}

bool Interface::isBoundaryInterface()
{
    return this->BoundaryConditionPointer != NULL;
}

string Interface::toString()
{
	ostringstream tmp;
    if (this->isBoundaryInterface())
    {
        tmp << "Boundary Interface at: \n";
        tmp << "(" << this->center.x << "," << this->center.y << ") \n";
        tmp << "Connected Cell: \n";
        if (negCell == NULL) 
            tmp << this->posCell->toString() << "\n";
        else
            tmp << this->negCell->toString() << "\n";
        tmp << "\n";
    }
    else
    {
        tmp << "Interface at: \n";
        tmp << "(" << this->center.x << "," << this->center.y << ") \n";
        tmp << "Connected Cell: \n";
        tmp << this->negCell->toString() << "\n";
        tmp << this->posCell->toString() << "\n";
        if (this->isGhostInterface())
            tmp << "GhostInterface\n";
        tmp << "\n";
    }
	return tmp.str();
}

string Interface::writeCenter()
{
    ostringstream tmp;
    tmp << this->center.x << " " << this->center.y << " 0.0\n";
    return tmp.str();
}

void Interface::interpolatePrim(double * prim)
{
    // This method computes the Values of the primary variables at the interface
    // with linear interpolation

    prim[0] = 0.5*( this->negCell->getPrim().rho
                  + this->posCell->getPrim().rho );
    prim[1] = 0.5*( this->negCell->getPrim().U
                  + this->posCell->getPrim().U   );
    prim[2] = 0.5*( this->negCell->getPrim().V
                  + this->posCell->getPrim().V   );
    prim[3] = 0.5*( this->negCell->getPrim().L
                  + this->posCell->getPrim().L   );
}

void Interface::interpolatePrimThirdOrder(double * prim)
{
    // For Boundary Interfaces use only linear Interpolation
    if ( this->negCell->isGhostCell() || this->posCell->isGhostCell() )
    {
        this->interpolatePrim(prim);
        return;
    }

    prim[0] = 7.0 / 12.0 * ( this->posCell->getPrim().rho                        + this->negCell->getPrim().rho )
            - 1.0 / 12.0 * ( this->posCell->getOpposingCell(this)->getPrim().rho + this->negCell->getOpposingCell(this)->getPrim().rho );

    prim[1] = 7.0 / 12.0 * ( this->posCell->getPrim().U                        + this->negCell->getPrim().U )
            - 1.0 / 12.0 * ( this->posCell->getOpposingCell(this)->getPrim().U + this->negCell->getOpposingCell(this)->getPrim().U );

    prim[2] = 7.0 / 12.0 * ( this->posCell->getPrim().V                        + this->negCell->getPrim().V )
            - 1.0 / 12.0 * ( this->posCell->getOpposingCell(this)->getPrim().V + this->negCell->getOpposingCell(this)->getPrim().V );

    prim[3] = 7.0 / 12.0 * ( this->posCell->getPrim().L                        + this->negCell->getPrim().L )
            - 1.0 / 12.0 * ( this->posCell->getOpposingCell(this)->getPrim().L + this->negCell->getOpposingCell(this)->getPrim().L );
}

void Interface::interpolateCons(double * cons)
{
    // This method computes the rhoValues of the consary variables at the interface
    // with linear interpolation

    cons[0] = 0.5*( this->negCell->getCons().rho
                  + this->posCell->getCons().rho );

    cons[1] = 0.5*( this->negCell->getCons().rhoU
                  + this->posCell->getCons().rhoU );

    cons[2] = 0.5*( this->negCell->getCons().rhoV
                  + this->posCell->getCons().rhoV );

    cons[3] = 0.5*( this->negCell->getCons().rhoE
                  + this->posCell->getCons().rhoE );
}

void Interface::interpolateConsThirdOrder(double * cons)
{
    // For Boundary Interfaces use only linear Interpolation
    if ( this->negCell->isGhostCell() || this->posCell->isGhostCell() )
    {
        this->interpolateCons(cons);
        return;
    }

    cons[0] = 7.0 / 12.0 * ( this->posCell->getCons().rho                        + this->negCell->getCons().rho )
            - 1.0 / 12.0 * ( this->posCell->getOpposingCell(this)->getCons().rho + this->negCell->getOpposingCell(this)->getCons().rho );

    cons[1] = 7.0 / 12.0 * ( this->posCell->getCons().rhoU                        + this->negCell->getCons().rhoU )
            - 1.0 / 12.0 * ( this->posCell->getOpposingCell(this)->getCons().rhoU + this->negCell->getOpposingCell(this)->getCons().rhoU );

    cons[2] = 7.0 / 12.0 * ( this->posCell->getCons().rhoV                        + this->negCell->getCons().rhoV )
            - 1.0 / 12.0 * ( this->posCell->getOpposingCell(this)->getCons().rhoV + this->negCell->getOpposingCell(this)->getCons().rhoV );

    cons[3] = 7.0 / 12.0 * ( this->posCell->getCons().rhoE                        + this->negCell->getCons().rhoE )
            - 1.0 / 12.0 * ( this->posCell->getOpposingCell(this)->getCons().rhoE + this->negCell->getOpposingCell(this)->getCons().rhoE );
}

void Interface::differentiateConsNormal(double* normalGradCons, double* prim)
{
    // This method computes the spacial derivatives of the conservative Variables.
    // The derivatives are computed by central finite differences.

    // ========================================================================
    // normal direction
    // ========================================================================

    // compute the distance between 
    double dn = ( ( this->posCell->getDx().x + this->negCell->getDx().x ) * normal.x
        + ( this->posCell->getDx().y + this->negCell->getDx().y ) * normal.y ) * 0.5;

    normalGradCons[0] = ( this->posCell->getCons().rho  - this->negCell->getCons().rho )  / (dn * prim[0] );

    normalGradCons[1] = ( this->posCell->getCons().rhoU - this->negCell->getCons().rhoU ) / ( dn * prim[0] );

    normalGradCons[2] = ( this->posCell->getCons().rhoV - this->negCell->getCons().rhoV ) / ( dn * prim[0] );

    normalGradCons[3] = ( this->posCell->getCons().rhoE - this->negCell->getCons().rhoE ) / ( dn * prim[0] );
}

void Interface::differentiateConsNormalThirdOrder(double* normalGradCons, double* prim)
{
    // This method computes the spacial derivatives of the conservative Variables.
    // The derivatives are computed by third order central finite differences.

    if ( this->negCell->isGhostCell() || this->posCell->isGhostCell() )
    {
        this->differentiateConsNormal(normalGradCons, prim);
        return;
    }

    // ========================================================================
    // normal direction
    // ========================================================================

    // compute the distance between 
    double dn = ( ( this->posCell->getDx().x + this->negCell->getDx().x ) * normal.x
        + ( this->posCell->getDx().y + this->negCell->getDx().y ) * normal.y ) * 0.5;

    normalGradCons[0] = ( 5.0/4.0  * ( this->posCell->getCons().rho                        - this->negCell->getCons().rho )  
                        - 1.0/12.0 * ( this->posCell->getOpposingCell(this)->getCons().rho - this->negCell->getOpposingCell(this)->getCons().rho )
                        ) / ( dn * prim[0] );

    normalGradCons[1] = ( 5.0/4.0  * ( this->posCell->getCons().rhoU                        - this->negCell->getCons().rhoU )  
                        - 1.0/12.0 * ( this->posCell->getOpposingCell(this)->getCons().rhoU - this->negCell->getOpposingCell(this)->getCons().rhoU )
                        ) / ( dn * prim[0] );

    normalGradCons[2] = ( 5.0/4.0  * ( this->posCell->getCons().rhoV                        - this->negCell->getCons().rhoV )  
                        - 1.0/12.0 * ( this->posCell->getOpposingCell(this)->getCons().rhoV - this->negCell->getOpposingCell(this)->getCons().rhoV )
                        ) / ( dn * prim[0] );

    normalGradCons[3] = ( 5.0/4.0  * ( this->posCell->getCons().rhoE                        - this->negCell->getCons().rhoE )  
                        - 1.0/12.0 * ( this->posCell->getOpposingCell(this)->getCons().rhoE - this->negCell->getOpposingCell(this)->getCons().rhoE )
                        ) / ( dn * prim[0] );

}

void Interface::differentiateConsTangential(double* tangentialGradCons, double* prim)
{
    // ========================================================================
    // tangential direction
    // ========================================================================

    // The tangential derivative is computed by finite difference between the
    // values at the edge of the interface (A, B in fig).
    // These are computed by interpolation.
    //
    //  A = 0.5 (pos + pos pos + neg + neg pos)
    //
    //  ---------------------------------
    //  |               |               |
    //  |    neg pos    |    pos pos    |
    //  |               |               |
    //  --------------- A --------------
    //  |               |               |
    //  |               |               |
    //  |      neg      |      pos      |
    //  |               |               |
    //  |               |               |
    //  --------------- B ---------------
    //  |               |               |
    //  |    neg neg    |    pos neg    |
    //  |               |               |
    //  ---------------------------------

    // get the indieces of the perpendicular interfaces for tangential derivative
    int posIdx;
    int negIdx;
    if (this->axis == 0)
    {
        posIdx = 3;
        negIdx = 1;
    }
    else
    {
        posIdx = 2;
        negIdx = 0;
    }

    // compute the tangential distance (length of the interface)
    double dt = this->posCell->getDx().x * normal.y
              + this->posCell->getDx().y * normal.x;

    tangentialGradCons[0] = ( ( this->posCell->getNeighborCell(posIdx)->getCons().rho
                              + this->negCell->getNeighborCell(posIdx)->getCons().rho
                              + this->posCell->getCons().rho 
                              + this->negCell->getCons().rho 
                              ) * 0.25
                            - ( this->posCell->getNeighborCell(negIdx)->getCons().rho
                              + this->negCell->getNeighborCell(negIdx)->getCons().rho 
                              + this->posCell->getCons().rho 
                              + this->negCell->getCons().rho
                              ) * 0.25
                            ) / ( dt * prim[0] );

    tangentialGradCons[1] = ( ( this->posCell->getNeighborCell(posIdx)->getCons().rhoU
                              + this->negCell->getNeighborCell(posIdx)->getCons().rhoU
                              + this->posCell->getCons().rhoU 
                              + this->negCell->getCons().rhoU 
                              ) * 0.25
                            - ( this->posCell->getNeighborCell(negIdx)->getCons().rhoU
                              + this->negCell->getNeighborCell(negIdx)->getCons().rhoU 
                              + this->posCell->getCons().rhoU 
                              + this->negCell->getCons().rhoU 
                              ) * 0.25
                            ) / ( dt * prim[0] );

    tangentialGradCons[2] = ( ( this->posCell->getNeighborCell(posIdx)->getCons().rhoV
                              + this->negCell->getNeighborCell(posIdx)->getCons().rhoV
                              + this->posCell->getCons().rhoV 
                              + this->negCell->getCons().rhoV 
                              ) * 0.25
                            - ( this->posCell->getNeighborCell(negIdx)->getCons().rhoV
                              + this->negCell->getNeighborCell(negIdx)->getCons().rhoV 
                              + this->posCell->getCons().rhoV 
                              + this->negCell->getCons().rhoV 
                              ) * 0.25
                            ) / ( dt * prim[0] );

    tangentialGradCons[3] = ( ( this->posCell->getNeighborCell(posIdx)->getCons().rhoE
                              + this->negCell->getNeighborCell(posIdx)->getCons().rhoE
                              + this->posCell->getCons().rhoE 
                              + this->negCell->getCons().rhoE 
                              ) * 0.25
                            - ( this->posCell->getNeighborCell(negIdx)->getCons().rhoE
                              + this->negCell->getNeighborCell(negIdx)->getCons().rhoE 
                              + this->posCell->getCons().rhoE 
                              + this->negCell->getCons().rhoE 
                              ) * 0.25
                            ) / ( dt * prim[0] );

}

void Interface::computeTimeDerivative(double * prim, double * MomentU, double * MomentV, double * MomentXi,
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

void Interface::assembleFlux(double * MomentU, double * MomentV, double * MomentXi, double * a, double * b, double * A, double * timeCoefficients, double dy, double* prim, double tau)
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

void Interface::rotate(double * vector)
{
    double tmp = vector[1];
    vector[1] = vector[2];
    vector[2] = tmp;
}

void Interface::cons2prim(double * prim, double * cons)
{
    prim[0] = cons[0];
    prim[1] = cons[1] / cons[0];
    prim[2] = cons[2] / cons[0];
    prim[3] = ( fluidParam.K + 2.0 )*cons[0]
        / ( 4.0 * (cons[3] - 0.5*( cons[1] * cons[1] + cons[2] * cons[2] ) / cons[0] ) );
}

double Interface::distance(float2 point)
{
    return sqrt( ( this->center.x - point.x )*( this->center.x - point.x )
               + ( this->center.y - point.y )*( this->center.y - point.y ) );
}

void Interface::computeMicroSlope(double * prim, double * macroSlope, double * microSlope)
{
    // this method computes the micro slopes from the slopes of the conservative variables
    // the resulting microslopes contain the density, since they are computed from the slopes
    // of the conservative variables, which are rho, rhoU, rhoV and rhoE

    double A, B, C, E;

    // ========================================================================
    // this is the total energy density E = rhoE / rho
    E = prim[1] * prim[1] + prim[2] * prim[2] + ( this->fluidParam.K + 2.0 ) / ( 2.0*prim[3] );
    // ========================================================================

    // ========================================================================
    // the product rule of derivations is used here!
    A = 2.0*macroSlope[3] - E       * macroSlope[0];    // = 2 rho dE/dx
    B =     macroSlope[1] - prim[1] * macroSlope[0];    // =   rho dU/dx
    C =     macroSlope[2] - prim[2] * macroSlope[0];    // =   rho dV/dx
    // ========================================================================

    // compute micro slopes of primitive variables from macro slopes of conservatice variables
    microSlope[3] = (4.0 * prim[3]*prim[3])/(this->fluidParam.K + 2.0)
                  * ( A - 2.0*prim[1]*B - 2.0*prim[2]*C );

    microSlope[2] = 2.0 * prim[3] * C - prim[2] * microSlope[3];

    microSlope[1] = 2.0 * prim[3] * B - prim[1] * microSlope[3];

    microSlope[0] = macroSlope[0] - prim[1]*microSlope[1] - prim[2]*microSlope[2] - 0.5 * E* microSlope[3];
}

void Interface::computeMoments(double * prim, double * MomentU, double* MomentV, double * MomentXi, int numberMoments)
{
    //==================== U Moments ==========================================
    MomentU[0] = 1.0;
    MomentU[1] = prim[1];
    for (int i = 2; i < numberMoments; i++)
        MomentU[i] = prim[1] * MomentU[i - 1] + (i - 1)/(2.0*prim[3])*MomentU[i - 2];

    //==================== V Moments ==========================================
    MomentV[0] = 1.0;
    MomentV[1] = prim[2];
    for (int i = 2; i < numberMoments; i++)
        MomentV[i] = prim[2] * MomentV[i - 1] + (i - 1)/(2.0*prim[3])*MomentV[i - 2];

    //==================== Xi Moments =========================================
    MomentXi[0] = 1.0;
    MomentXi[1] = 0.0;
    MomentXi[2] = this->fluidParam.K / (2.0 * prim[3]);
    MomentXi[3] = 0.0;
    MomentXi[4] = ( 2.0*this->fluidParam.K + 1.0*this->fluidParam.K*this->fluidParam.K ) / (4.0 * prim[3] * prim[3]);
    MomentXi[5] = 0.0;
    MomentXi[6] = 0.0;
}

