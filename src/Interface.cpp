

#include "Interface.h"
#include "IncompressibleInterface.h"
#include "CompressibleInterface.h"
#include "Types.h"
#include <sstream>

int Interface::interpolationOrder = 3;

Interface::Interface()
{
}

Interface::Interface(Cell* negCell, Cell* posCell, bool negAdd, bool posAdd, float2** nodes, FluidParameter fluidParam, BoundaryCondition* BC, double periodicLength)
{
    // ========================================================================
    //                  Copy attributes
    // ========================================================================
	this->negCell = negCell;
	this->posCell = posCell;

    this->nodes[0] = nodes[0];
    this->nodes[1] = nodes[1]; 

    this->BoundaryConditionPointer = BC;

    this->fluidParam = fluidParam;
    // ========================================================================
    
    // ========================================================================
    //                  Compute Center
    // ========================================================================
    this->center.x = 0.5 * (this->nodes[0]->x + this->nodes[1]->x);
    this->center.y = 0.5 * (this->nodes[0]->y + this->nodes[1]->y);
    // ========================================================================
    
    // ========================================================================
    //                  Compute interface area
    // ========================================================================
    this->area = sqrt( (this->nodes[0]->x - this->nodes[1]->x )*(this->nodes[0]->x - this->nodes[1]->x )
                     + (this->nodes[0]->y - this->nodes[1]->y )*(this->nodes[0]->y - this->nodes[1]->y ));
    // ========================================================================
    
    // ========================================================================
    //                  Compute Normal
    // ========================================================================
    //      -----[1]-------
    //            |             The normal is computed such that it points
    //            |   n         to the right when the one follows the
    //            |----->       vector  from the first to the second node.
    //    negCell | posCell     
    //            |             n = (N1 - N0) x (0 0 1)
    //      -----[0]-------    
    // ========================================================================        
    this->normal.x = - ( this->nodes[1]->y - this->nodes[0]->y ) / area;
    this->normal.y =   ( this->nodes[1]->x - this->nodes[0]->x ) / area;
    // ======================================================================== 
    
    // ========================================================================
    //                  Introduce Interfaces to Cells
    // ========================================================================
	if(this->negCell != NULL && negAdd) this->negCell->addInterface(this);
	if(this->posCell != NULL && posAdd) this->posCell->addInterface(this);
    // ========================================================================
    
    // ========================================================================
    //                  compute distances
    // ========================================================================
    // Ghostcells have the same distance as their pendants in the domain
    if(negCell != NULL){
        if(negAdd) this->negDistance =                  this->distance( this->negCell->getCenter() );
        else       this->negDistance = periodicLength - this->distance( this->negCell->getCenter() );
    }else{
        this->negDistance = this->distance( this->posCell->getCenter() );
    }
    if(posCell != NULL){
        if(posAdd) this->posDistance =                  this->distance( this->posCell->getCenter() );
        else       this->posDistance = periodicLength - this->distance( this->posCell->getCenter() );
    }else{
        this->posDistance = this->distance( this->negCell->getCenter() );
    }
    // ========================================================================
    
    // ========================================================================
    //                  Initialize Fluxes
    // ========================================================================
    for ( int i = 0; i < 4; i++ )
    {
        this->timeIntegratedFlux[i] = 0.0;
        this->FluxDensity[i] = 0.0;
    }
    // ========================================================================
}

Interface::~Interface()
{
}

// ============================================================================
//                     Interface Factory
// ============================================================================
Interface * Interface::createInterface(InterfaceType type, Cell * negCell, Cell * posCell, bool negAdd, bool posAdd, 
                                       float2** nodes, FluidParameter fluidParam, BoundaryCondition * BC, double periodicLength)
{
    Interface* tmp = NULL;

    if ( type == incompressible )
        tmp = new IncompressibleInterface(negCell, posCell, negAdd, posAdd, nodes, fluidParam, BC, periodicLength);
    else
        tmp = new CompressibleInterface(negCell, posCell, negAdd, posAdd, nodes, fluidParam, BC, periodicLength);

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
    double normalGradConsTest[4];
    double tangentialGradCons[4];
    double timeGrad[4];

    double a[4];
    double b[4] = {0.0, 0.0, 0.0, 0.0};
    double A[4];

    double MomentU[NUMBER_OF_MOMENTS];
    double MomentV[NUMBER_OF_MOMENTS];
    double MomentXi[NUMBER_OF_MOMENTS];

    // ========================================================================
    // interpolated primary variables at the interface
    this->interpolatePrim(prim);

    // spacial gradients of the conservative varibles
    //this->differentiateConsNormal(normalGradCons, prim);
    this->differentiateConsLeastSquare(normalGradCons, tangentialGradCons, prim);
    // ========================================================================
    
    // ========================================================================
    //          Some Tests
    // ========================================================================
    PrimitiveVariable posPrim = this->posCell->getPrim();
    PrimitiveVariable negPrim = this->negCell->getPrim();
    this->transformGlobal2Local( (double*)&posPrim );
    this->transformGlobal2Local( (double*)&negPrim );
    // ========================================================================
    
    // ========================================================================
    // Transformation in local coordinate system
    // ========================================================================
    transformGlobal2Local(prim);
    transformGlobal2Local(normalGradCons);
    // ========================================================================
    
    // ========================================================================
    // spacial micro slopes a = a1 + a2 u + a3 v + 0.5 a4 (u^2 + v^2 + xi^2)
    // ========================================================================
    this->computeMicroSlope(prim, normalGradCons, a);
    this->computeMicroSlope(prim, tangentialGradCons, b);
    // ========================================================================

    // ========================================================================
    this->computeMoments(prim, MomentU, MomentV, MomentXi, NUMBER_OF_MOMENTS);
    // ========================================================================

    // ========================================================================
    // temporal micro slopes A = A1 + A2 u + A3 v + 0.5 A4 (u^2 + v^2 + xi^2)
    // ========================================================================
    this->computeTimeDerivative(prim, MomentU, MomentV, MomentXi, a, b, timeGrad);

    this->computeMicroSlope(prim, timeGrad, A);
    // ========================================================================

    // ========================================================================
    // Formular as in the Rayleigh-Bernard-Paper (Xu, Lui, 1999)
    // ========================================================================
    double tau = 2.0*prim[3] * this->fluidParam.nu;
    // ========================================================================
    
    // ========================================================================
    // time integration Coefficients
    // ========================================================================
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };
    // ========================================================================

    // ========================================================================
    // compute mass and momentum fluxes
    // ========================================================================
    this->assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients, prim, tau);
    // ========================================================================
    
    transformLocal2Global(this->timeIntegratedFlux);
    transformLocal2Global(this->timeIntegratedFlux_1);
    transformLocal2Global(this->timeIntegratedFlux_2);
    transformLocal2Global(this->timeIntegratedFlux_3);
    transformLocal2Global(this->FluxDensity);
    
    // ========================================================================
    // add Flux values to cell updates
    // ========================================================================
    //this->posCell->addFlux(this->timeIntegratedFlux,  1.0, this);
    //this->negCell->addFlux(this->timeIntegratedFlux, -1.0, this);
    // ========================================================================

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

ConservedVariable Interface::getTimeIntegratedFlux_1()
{
    ConservedVariable tmp;
    tmp.rho  = this->timeIntegratedFlux_1[0];
    tmp.rhoU = this->timeIntegratedFlux_1[1];
    tmp.rhoV = this->timeIntegratedFlux_1[2];
    tmp.rhoE = this->timeIntegratedFlux_1[3];
    return tmp;
}

ConservedVariable Interface::getTimeIntegratedFlux_2()
{
    ConservedVariable tmp;
    tmp.rho  = this->timeIntegratedFlux_2[0];
    tmp.rhoU = this->timeIntegratedFlux_2[1];
    tmp.rhoV = this->timeIntegratedFlux_2[2];
    tmp.rhoE = this->timeIntegratedFlux_2[3];
    return tmp;
}

ConservedVariable Interface::getTimeIntegratedFlux_3()
{
    ConservedVariable tmp;
    tmp.rho  = this->timeIntegratedFlux_3[0];
    tmp.rhoU = this->timeIntegratedFlux_3[1];
    tmp.rhoV = this->timeIntegratedFlux_3[2];
    tmp.rhoE = this->timeIntegratedFlux_3[3];
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
    double sign = ( askingCell->getCenter().x - this->center.x ) * this->normal.x
                + ( askingCell->getCenter().y - this->center.y ) * this->normal.y ;

    return sign/fabs(sign);
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

float2 Interface::getNormal()
{
    return normal;
}

float2 Interface::getCenter()
{
    return this->center;
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

float2 Interface::getPosConnectivity()
{
    float2 vector;
    vector.x = this->posCell->getCenter().x - this->center.x;
    vector.y = this->posCell->getCenter().y - this->center.y;
    return vector;
}

float2 Interface::getNegConnectivity()
{
    float2 vector;
    vector.x = this->negCell->getCenter().x - this->center.x;
    vector.y = this->negCell->getCenter().y - this->center.y;
    return vector;
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

    double distance = this->posDistance + this->negDistance;

    prim[0] = ( this->negCell->getPrim().rho * this->posDistance
              + this->posCell->getPrim().rho * this->negDistance )
            / ( distance );

    prim[1] = ( this->negCell->getPrim().U   * this->posDistance
              + this->posCell->getPrim().U   * this->negDistance )
            / ( distance );

    prim[2] = ( this->negCell->getPrim().V   * this->posDistance
              + this->posCell->getPrim().V   * this->negDistance )
            / ( distance );

    prim[3] = ( this->negCell->getPrim().L   * this->posDistance
              + this->posCell->getPrim().L   * this->negDistance )
            / ( distance );

    int i = 0;
}

void Interface::differentiateConsNormal(double* normalGradCons, double* prim)
{
    // This method computes the spacial derivatives of the conservative Variables.
    // The derivatives are computed by central finite differences.

    // ========================================================================
    // normal direction
    // ========================================================================

    // compute the distance between the adjacent cell centers
    double distance = this->posDistance + this->negDistance;

    normalGradCons[0] = ( this->posCell->getCons().rho  - this->negCell->getCons().rho )  / ( distance * prim[0] );

    normalGradCons[1] = ( this->posCell->getCons().rhoU - this->negCell->getCons().rhoU ) / ( distance * prim[0] );

    normalGradCons[2] = ( this->posCell->getCons().rhoV - this->negCell->getCons().rhoV ) / ( distance * prim[0] );

    normalGradCons[3] = ( this->posCell->getCons().rhoE - this->negCell->getCons().rhoE ) / ( distance * prim[0] );
}

void Interface::differentiateConsLeastSquare(double* normalGradCons, double* tangentialGradCons, double* prim)
{
    // ========================================================================
    // Read gradients from adjacent Cells and transform into local coordinates
    // ========================================================================

    ConservedVariable posGradN;
    posGradN.rho  =   this->normal.x * this->posCell->getGradientX().rho  + this->normal.y * this->posCell->getGradientY().rho ;
    posGradN.rhoU =   this->normal.x * this->posCell->getGradientX().rhoU + this->normal.y * this->posCell->getGradientY().rhoU;
    posGradN.rhoV =   this->normal.x * this->posCell->getGradientX().rhoV + this->normal.y * this->posCell->getGradientY().rhoV;
    posGradN.rhoE =   this->normal.x * this->posCell->getGradientX().rhoE + this->normal.y * this->posCell->getGradientY().rhoE;
    ConservedVariable posGradT;
    posGradT.rho  = - this->normal.y * this->posCell->getGradientX().rho  + this->normal.x * this->posCell->getGradientY().rho ;
    posGradT.rhoU = - this->normal.y * this->posCell->getGradientX().rhoU + this->normal.x * this->posCell->getGradientY().rhoU;
    posGradT.rhoV = - this->normal.y * this->posCell->getGradientX().rhoV + this->normal.x * this->posCell->getGradientY().rhoV;
    posGradT.rhoE = - this->normal.y * this->posCell->getGradientX().rhoE + this->normal.x * this->posCell->getGradientY().rhoE;

    ConservedVariable negGradN;
    negGradN.rho  =   this->normal.x * this->negCell->getGradientX().rho  + this->normal.y * this->negCell->getGradientY().rho ;
    negGradN.rhoU =   this->normal.x * this->negCell->getGradientX().rhoU + this->normal.y * this->negCell->getGradientY().rhoU;
    negGradN.rhoV =   this->normal.x * this->negCell->getGradientX().rhoV + this->normal.y * this->negCell->getGradientY().rhoV;
    negGradN.rhoE =   this->normal.x * this->negCell->getGradientX().rhoE + this->normal.y * this->negCell->getGradientY().rhoE;
    ConservedVariable negGradT;
    negGradT.rho  = - this->normal.y * this->negCell->getGradientX().rho  + this->normal.x * this->negCell->getGradientY().rho ;
    negGradT.rhoU = - this->normal.y * this->negCell->getGradientX().rhoU + this->normal.x * this->negCell->getGradientY().rhoU;
    negGradT.rhoV = - this->normal.y * this->negCell->getGradientX().rhoV + this->normal.x * this->negCell->getGradientY().rhoV;
    negGradT.rhoE = - this->normal.y * this->negCell->getGradientX().rhoE + this->normal.x * this->negCell->getGradientY().rhoE;
    // ========================================================================
    
    // ========================================================================
    // Average Gradients
    // ========================================================================
    double distance = this->posDistance + this->negDistance;

    normalGradCons[0] = ( posGradN.rho  * this->posDistance + negGradN.rho  * this->negDistance ) / ( distance * prim[0] );
    normalGradCons[1] = ( posGradN.rhoU * this->posDistance + negGradN.rhoU * this->negDistance ) / ( distance * prim[0] );
    normalGradCons[2] = ( posGradN.rhoV * this->posDistance + negGradN.rhoV * this->negDistance ) / ( distance * prim[0] );
    normalGradCons[3] = ( posGradN.rhoE * this->posDistance + negGradN.rhoE * this->negDistance ) / ( distance * prim[0] );

    tangentialGradCons[0] = ( posGradT.rho  * this->posDistance + negGradT.rho  * this->negDistance ) / ( distance * prim[0] );
    tangentialGradCons[1] = ( posGradT.rhoU * this->posDistance + negGradT.rhoU * this->negDistance ) / ( distance * prim[0] );
    tangentialGradCons[2] = ( posGradT.rhoV * this->posDistance + negGradT.rhoV * this->negDistance ) / ( distance * prim[0] );
    tangentialGradCons[3] = ( posGradT.rhoE * this->posDistance + negGradT.rhoE * this->negDistance ) / ( distance * prim[0] );
    // ========================================================================

    int i = 0;
}

void Interface::transformGlobal2Local(double * vec)
{
    // euclidian components in global coordinatesystem
    double u0 = vec[1];
    double v0 = vec[2];

    // transformation in local coordinatesystem
    // n = (n1,n2)
    // t = (-n2,n1)
    // vL = [n t]^T * v0
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
    // vL = [n t] * v0
    vec[1] = this->normal.x * un - this->normal.y * vt;
    vec[2] = this->normal.y * un + this->normal.x * vt;
}

double Interface::distance(float2 point)
{
    // Compute the projected distance of a point on the normal of the interface
    return fabs( ( this->center.x - point.x )*( this->normal.x )
               + ( this->center.y - point.y )*( this->normal.y ) );
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

