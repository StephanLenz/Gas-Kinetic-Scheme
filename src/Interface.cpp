

#include "Interface.h"
#include "IncompressibleInterface.h"
#include "CompressibleInterface.h"
#include "Types.h"
#include <sstream>

int Interface::interpolationOrder = 3;

Interface::Interface()
{
}

Interface::Interface(Cell* negCell, Cell* posCell, float2** nodes, FluidParameter fluidParam, BoundaryCondition* BC)
{
	this->negCell = negCell;
	this->posCell = posCell;

    this->nodes[0] = nodes[0];
    this->nodes[1] = nodes[1];

    this->center.x = 0.5 * (this->nodes[0]->x + this->nodes[1]->x);
    this->center.y = 0.5 * (this->nodes[0]->y + this->nodes[1]->y);

    this->area = sqrt( (this->nodes[0]->x - this->nodes[1]->x )*(this->nodes[0]->x - this->nodes[1]->x )
                     + (this->nodes[0]->y - this->nodes[1]->y )*(this->nodes[0]->y - this->nodes[1]->y ));

    this->normal.x = - ( this->nodes[1]->y - this->nodes[0]->y ) / area;
    this->normal.y =   ( this->nodes[1]->x - this->nodes[0]->x ) / area;

    this->BoundaryConditionPointer = BC;
    
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

	if(this->negCell != NULL) this->negCell->addInterface(this);
	if(this->posCell != NULL) this->posCell->addInterface(this);

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
Interface * Interface::createInterface(InterfaceType type, Cell * negCell, Cell * posCell, float2** nodes, FluidParameter fluidParam, BoundaryCondition * BC)
{
    Interface* tmp = NULL;

    if ( type == incompressible )
        tmp = new IncompressibleInterface(negCell, posCell, nodes, fluidParam, BC);
    else
        tmp = new CompressibleInterface(negCell, posCell, nodes, fluidParam, BC);

    return tmp;
}

void Interface::computeFlux(double dt)
{
    this->computeInternalFlux(dt);
}

void Interface::computeInternalFlux(double dt)
{
    const int NUMBER_OF_MOMENTS = 7;

    double prim[4];
    double normalGradCons[4];
    double timeGrad[4];

    double a[4];
    double b[] = {0.0, 0.0, 0.0, 0.0};
    double A[4];

    double MomentU[NUMBER_OF_MOMENTS];
    double MomentV[NUMBER_OF_MOMENTS];
    double MomentXi[NUMBER_OF_MOMENTS];

    // ========================================================================
    // interpolated primary variables at the interface
    this->interpolatePrim(prim);

    // spacial gradients of the conservative varibles
    this->differentiateConsNormal(normalGradCons, prim);
    // ========================================================================

    
    // ========================================================================
    // Transformation in local coordinate system
    transformGlobal2Local(prim);
    transformGlobal2Local(normalGradCons);
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
    // spacial micro slopes a = a1 + a2 u + a3 v + 0.4 a4 (u^2 + v^2 + xi^2)
    this->computeMicroSlope(prim, normalGradCons, a);
    // ========================================================================

    // ========================================================================
    this->computeMoments(prim, MomentU, MomentV, MomentXi, NUMBER_OF_MOMENTS);
    // ========================================================================

    // ========================================================================
    // temporal micro slopes A = A1 + A2 u + A3 v + 0.4 A4 (u^2 + v^2 + xi^2)
    this->computeTimeDerivative(prim, MomentU, MomentV, MomentXi, a, b, timeGrad);

    this->computeMicroSlope(prim, timeGrad, A);
    // ========================================================================

    // ========================================================================
    // compute mass and momentum fluxes
    this->assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients, prim, tau);
    // ========================================================================
    
    transformLocal2Global(this->timeIntegratedFlux);
    transformGlobal2Local(this->FluxDensity);

    int i = 1;
}

Cell * Interface::getNeigborCell(Cell * askingCell)
{
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

double Interface::getFluxSign(Cell * askingCell)
{
    if      (askingCell == this->posCell) return  1.0;
    else if (askingCell == this->negCell) return -1.0;
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

void Interface::addCell(Cell * that)
{
    if(this->posCell == NULL)
        this->posCell = that;
    else if(this->negCell == NULL)
        this->negCell = that;
}

float2 * Interface::getNode(int i)
{
    return this->nodes[i];
}

float2 Interface::getScaledNormal()
{
    if(this->posCell == NULL)
        return float2(   this->normal.x * this->distance(negCell->getCenter()),   this->normal.y * this->distance(negCell->getCenter()) );
    else
        return float2( - this->normal.x * this->distance(posCell->getCenter()), - this->normal.y * this->distance(posCell->getCenter()) );
}

BoundaryCondition * Interface::getBoundaryCondition()
{
    return this->BoundaryConditionPointer;
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

    double negDistance = this->distance( this->negCell->getCenter() );
    double posDistance = this->distance( this->posCell->getCenter() );

    prim[0] = ( this->negCell->getPrim().rho * posDistance
              + this->posCell->getPrim().rho * negDistance)
            / ( negDistance + posDistance );

    prim[1] = ( this->negCell->getPrim().U   * posDistance
              + this->posCell->getPrim().U   * negDistance)
            / ( negDistance + posDistance );

    prim[2] = ( this->negCell->getPrim().V   * posDistance
              + this->posCell->getPrim().V   * negDistance)
            / ( negDistance + posDistance );

    prim[3] = ( this->negCell->getPrim().L   * posDistance
              + this->posCell->getPrim().L   * negDistance)
            / ( negDistance + posDistance );
}

void Interface::differentiateConsNormal(double* normalGradCons, double* prim)
{
    // This method computes the spacial derivatives of the conservative Variables.
    // The derivatives are computed by central finite differences.

    // ========================================================================
    // normal direction
    // ========================================================================

    // compute the distance between the adjacent cell centers

    double distance = this->distance( this->negCell->getCenter() )
                    + this->distance( this->posCell->getCenter() );

    normalGradCons[0] = ( this->posCell->getCons().rho  - this->negCell->getCons().rho )  / ( distance * prim[0] );

    normalGradCons[1] = ( this->posCell->getCons().rhoU - this->negCell->getCons().rhoU ) / ( distance * prim[0] );

    normalGradCons[2] = ( this->posCell->getCons().rhoV - this->negCell->getCons().rhoV ) / ( distance * prim[0] );

    normalGradCons[3] = ( this->posCell->getCons().rhoE - this->negCell->getCons().rhoE ) / ( distance * prim[0] );
}

void Interface::transformGlobal2Local(double * vec)
{
    // euclidian components in global coordinatesystem
    double u0 = vec[1];
    double v0 = vec[2];

    // transformation in local coordinatesystem
    // n = (n1,n2)
    // t = (-n2,n1)
    vec[1] =   this->normal.x * u0 + this->normal.y * v0;
    vec[2] = - this->normal.y * u0 + this->normal.x * v0;
}

void Interface::transformLocal2Global(double * vec)
{
    // euclidian components in local coordinatesystem
    double un = vec[1];
    double vt = vec[2];

    // transformation in global coordinatesystem
    // n = (n1,n2)
    // t = (-n2,n1)
    vec[1] = this->normal.x * un - this->normal.y * vt;
    vec[2] = this->normal.y * un + this->normal.x * vt;
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

