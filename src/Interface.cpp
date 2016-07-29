

#include "Interface.h"
#include "IncompressibleInterface.h"
#include "CompressibleInterface.h"
#include "Types.h"
#include <sstream>

int Interface::interpolationOrder = 3;

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

// ============================================================================
//                     Interface Factory
// ============================================================================
Interface * Interface::createInterface(InterfaceType type, Cell * negCell, Cell * posCell, float2 center, float2 normal, FluidParameter fluidParam, InterfaceBC * BC)
{
    Interface* tmp = NULL;

    if ( type == incompressible )
        tmp = new IncompressibleInterface(negCell, posCell, center, normal, fluidParam, BC);
    else
        tmp = new CompressibleInterface(negCell, posCell, center, normal, fluidParam, BC);

    return tmp;
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
    if( interpolationOrder == 1 )
    {
        this->interpolatePrim(prim);
    }
    else
    {
        this->interpolatePrimThirdOrder(prim);
    }

    // spacial gradients of the conservative varibles
    if( interpolationOrder == 1 )
        this->differentiateConsNormal(normalGradCons, prim);
    else
        this->differentiateConsNormalThirdOrder(normalGradCons, prim);

    this->differentiateConsTangential(tangentialGradCons, prim);
    // ========================================================================

    // ========================================================================
    // Formular as in the Rayleigh-Bernard-Paper (Xu, Lui, 1999)
    double tau = 2.0*prim[3] * this->fluidParam.nu;
    // ========================================================================
    
    // ========================================================================
    // time integration Coefficients
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };
    // ========================================================================
    
    // ========================================================================
    // spacial micro slopes a = a1 + a2 u + a3 v
    //                      b = b1 + b2 u + b3 v
    if ( this->axis == 0 )
    {
        this->computeMicroSlope(prim, normalGradCons, a);
        this->computeMicroSlope(prim, tangentialGradCons, b);
    }
    else
    {
        this->computeMicroSlope(prim, normalGradCons, b);
        this->computeMicroSlope(prim, tangentialGradCons, a);
    }
    // ========================================================================

    // This Block can turn off the usage of tangential Derivatives
    // ========================================================================
    // ========================================================================
    // ========================================================================
    for ( int i = 0; i < 4; i++ )
    {
        if ( this->axis == 0 )
            b[i] = 0.0;
        else
            a[i] = 0.0;
    }
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
    this->computeTimeDerivative(prim, MomentU, MomentV, MomentXi, a, b, timeGrad);

    this->computeMicroSlope(prim, timeGrad, A);
    // ========================================================================

    // ========================================================================
    // compute mass and momentum fluxes
    this->assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients, dy, prim, tau);
    // ========================================================================
    
    int i = 1;
}

void Interface::computeBoundaryFlux(double dt)
{
    const int NUMBER_OF_MOMENTS = 7;
    double timeGrad[4];

    double a[] = {0.0,0.0,0.0};
    double b[] = {0.0,0.0,0.0};
    double A[] = {0.0,0.0,0.0};

    double MomentU[NUMBER_OF_MOMENTS];
    double MomentV[NUMBER_OF_MOMENTS];
    double MomentXi[NUMBER_OF_MOMENTS];

    // ========================================================================
    double sign = 1.0;
    if ( posCell == NULL )
        sign = -1.0;

    // compute the length of the interface
    double dy = this->getCellInDomain()->getDx().x * normal.y
              + this->getCellInDomain()->getDx().y * normal.x;

    double distance = this->distance(this->getCellInDomain()->getCenter());
    // ========================================================================

    // ========================================================================
    double x1 = this->distance(this->getCellInDomain()->getCenter());
    double x2 = this->distance(this->getCellInDomain()->getOpposingCell(this)->getCenter());
    
    PrimaryVariable Z1 = this->getCellInDomain()->getPrim();
    PrimaryVariable Z2 = this->getCellInDomain()->getOpposingCell(this)->getPrim();
    // ========================================================================

    // ========================================================================
    PrimaryVariable prim;
    prim.rho = ( x2*x2*Z1.rho - x1*x1*Z2.rho )/( x2*x2 - x1*x1 );
    prim.U   = 0.0;
    prim.V   = 0.0;
    prim.L   = ( x2*x2*Z1.L - x1*x1*Z2.L )/( x2*x2 - x1*x1 );

    ConservedVariable normalGradCons;
    normalGradCons.rho  = 0.0;
    normalGradCons.rhoU = sign * ( x2*x2*Z1.U - x1*x1*Z2.U )/( (x2*x2 - x1*x1)*x1*x2 );
    normalGradCons.rhoV = sign * ( x2*x2*Z1.V - x1*x1*Z2.V )/( (x2*x2 - x1*x1)*x1*x2 );
    normalGradCons.rhoE = 0.0;
    // ========================================================================

    //// ========================================================================
    //PrimaryVariable prim;
    //prim.rho = this->getCellInDomain()->getPrim().rho;
    //prim.U   = 0.0;
    //prim.V   = 0.0;
    //prim.L   = this->getCellInDomain()->getPrim().L;

    //ConservedVariable normalGradCons;
    //normalGradCons.rho  = 0.0;
    //normalGradCons.rhoU = sign * this->getCellInDomain()->getCons().rhoU / (distance * prim.rho);
    //normalGradCons.rhoV = sign * this->getCellInDomain()->getCons().rhoV / (distance * prim.rho);
    //normalGradCons.rhoE = 0.0;
    //// ========================================================================

    // ========================================================================
    // Formular as in the Rayleigh-Bernard-Paper (Xu, Lui, 1999)
    double tau = 2.0 * prim.L * this->fluidParam.nu;
    // ========================================================================
    
    // ========================================================================
    // time integration Coefficients
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };
    // ========================================================================
 
    // ========================================================================
    // spacial micro slopes a = a1 + a2 u + a3 v
    //                      b = b1 + b2 u + b3 v
    if ( this->axis == 0 )
    {
        this->computeMicroSlope((double*)&prim, (double*)&normalGradCons, a);
    }
    else
    {
        this->computeMicroSlope((double*)&prim, (double*)&normalGradCons, b);
    }
    // ========================================================================
    
    // ========================================================================
    // The Moments are computed without the density
    // Therefore their dimensions are powers of velocity
    this->computeMoments((double*)&prim, MomentU, MomentV, MomentXi, NUMBER_OF_MOMENTS);
    // ========================================================================

    // ========================================================================
    // temporal micro slopes A = A1 + A2 u + A3 v
    // The temporal macro slopes (timeGrad) also contain the density by explicit multiplication
    this->computeTimeDerivative((double*)&prim, MomentU, MomentV, MomentXi, a, b, timeGrad);

    this->computeMicroSlope((double*)&prim, timeGrad, A);
    // ========================================================================

    // ========================================================================
    // compute mass and momentum fluxes
    this->assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients, dy, (double*)&prim, tau);
    // ========================================================================

    int i = 1;
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
    // ========================================================================
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

void Interface::setInterpolationOrder(int arg)
{
    interpolationOrder = arg;
}

int Interface::getInterpolationOrder()
{
    return interpolationOrder;
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
    // For Boundary Interfaces use only linear Interpolation (in case of GhostCells)
    if ( this->negCell->isGhostCell() || this->posCell->isGhostCell() )
    {
        this->interpolatePrim(prim);
        return;
    }

    // For first Interfaces in Domain use only linear Interpolation (in case of InterfaceBCs)
    if ( this->negCell->getOpposingCell(this) == NULL || this->posCell->getOpposingCell(this) == NULL )
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

    normalGradCons[0] = ( this->posCell->getCons().rho  - this->negCell->getCons().rho )  / ( dn * prim[0] );

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

void Interface::computeMoments(double * prim, double * MomentU, double* MomentV, double * MomentXi, int numberMoments)
{
    //==================== U Moments ==========================================
    MomentU[0] = 1.0;
    MomentU[1] = prim[1];
    for (int i = 2; i < numberMoments; i++)
        MomentU[i] = prim[1] * MomentU[i - 1] + ((i - 1)*MomentU[i - 2])/(2.0*prim[3]);

    //==================== V Moments ==========================================
    MomentV[0] = 1.0;
    MomentV[1] = prim[2];
    for (int i = 2; i < numberMoments; i++)
        MomentV[i] = prim[2] * MomentV[i - 1] + ((i - 1)*MomentV[i - 2])/(2.0*prim[3]);

    //==================== Xi Moments =========================================
    MomentXi[0] = 1.0;
    MomentXi[1] = 0.0;
    MomentXi[2] = this->fluidParam.K / (2.0 * prim[3]);
    MomentXi[3] = 0.0;
    MomentXi[4] = ( 2.0*this->fluidParam.K + 1.0*this->fluidParam.K*this->fluidParam.K ) / (4.0 * prim[3] * prim[3]);
    MomentXi[5] = 0.0;
    MomentXi[6] = ( 1.0*this->fluidParam.K + 4.0 ) / ( 2.0 * prim[3] ) * MomentXi[4];
}

