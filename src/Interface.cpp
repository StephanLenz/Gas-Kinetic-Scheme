// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      Interface.h
//
//          This class in an abstract baseclass for compressible and
//          incompressible Interfaces. It implements all methods that do not
//          depend on the type of interface.
//
//      Function:
//          Storage of Interface data
//          Implementation of Interface related computations
//              Flux computation 
//
// ============================================================================

#include "Interface.h"
#include "IncompressibleInterface.h"
#include "CompressibleInterface.h"
#include "Types.h"
#include <sstream>

int Interface::interpolationOrder = 1;
unsigned long int Interface::numberOfCells = 1;

Interface::Interface()
{
}

// ============================================================================
//      This constructor initializes the Interface.
//      It also computes:
//          Center
//          Area
//          Normal
//          Distance to pos and neg Cell
//
//      Parameters:
//          negCell:            Pointer to the adjacent cell in neg direction
//          posCell:            Pointer to the adjacent cell in pos direction
//          negAdd:             true if this interface should be introduced to 
//                              the neg Cell
//          posAdd:             true if this interface should be introduced to
//                              the pos Cell
//          nodes:              pointer to an array containing pointers to the
//                              two nodes defining this interface
//          fluidParam:         FluidParameter Object
//          BC:                 Pointer to Boundary Condition Object
//          periodicLengthX:    length of the domain in X direction
//          periodicLengthY:    length of the domain in Y direction
// ============================================================================
Interface::Interface(Cell* negCell, Cell* posCell, bool negAdd, bool posAdd, Node** nodes, FluidParameter fluidParam, BoundaryCondition* BC, double periodicLengthX, double periodicLengthY)
{
    this->ID = Interface::numberOfCells++;

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
    //                  compute distances from Interface center to Cell center
    // ========================================================================
    // This Method is called before the creation of the ghostcells. Ghost cells
    // have the same distance as their pendants in the domain.
    //
    // In the case of periodic boundaries, the pure distance would not be
    // correct. Therefore if the interface has a neighbor cell that should
    // not be added (because it lies on the opposite side of the domain), the 
    // distance is computed as the complemantary distance.
    // ========================================================================
    if(negCell != NULL){
        if(negAdd) this->negDistance =                   this->distance( this->negCell->getCenter() );
        else       this->negDistance = periodicLengthX - this->distance( this->negCell->getCenter() );
    }else{
        this->negDistance = this->distance( this->posCell->getCenter() );
    }
    if(posCell != NULL){
        if(posAdd) this->posDistance =                   this->distance( this->posCell->getCenter() );
        else       this->posDistance = periodicLengthX - this->distance( this->posCell->getCenter() );
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
//
//          This method should be used to create polymorphic Interfaces.
//
// ============================================================================
Interface * Interface::createInterface(InterfaceType type, Cell * negCell, Cell * posCell, bool negAdd, bool posAdd, 
                                       Node** nodes, FluidParameter fluidParam, BoundaryCondition * BC, double periodicLengthX, double periodicLengthY)
{
    Interface* tmp = NULL;

    if ( type == incompressible )
        tmp = new IncompressibleInterface(negCell, posCell, negAdd, posAdd, nodes, fluidParam, BC, periodicLengthX, periodicLengthY);
    else
        tmp = new CompressibleInterface(negCell, posCell, negAdd, posAdd, nodes, fluidParam, BC, periodicLengthX, periodicLengthY);

    return tmp;
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          Computation Methods
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method is currently not used
// ============================================================================
void Interface::setInterpolationOrder(int arg)
{
    interpolationOrder = arg;
}

// ============================================================================
//      This method is currently not used
// ============================================================================
int Interface::getInterpolationOrder()
{
    return interpolationOrder;
}

// ============================================================================
//      This method adds a connectivity from an interface to a cell on the 
//      side of the interface that is not yet occupied.
// ============================================================================
void Interface::addCell(Cell * newCell)
{
    if(this->posCell == NULL)
    {
        this->posCell = newCell;
        this->posDistance = this->distance( this->posCell->getCenter() );
    }
    else if(this->negCell == NULL)
    {
        this->negCell = newCell;
        this->negDistance = this->distance( this->negCell->getCenter() );
    }
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          Computation Methods
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method is the public interface to start the computation of the
//      Fluxes.
//
//      Parameters:
//          dt:   time step
// ============================================================================
void Interface::computeFlux(double dt)
{
    this->computeInternalFlux(dt);
}

// ====================================================================================================================
//
//                          Reconstruction of primitive variables
//
// ====================================================================================================================

// ============================================================================
//      This method computes the Values of the primary variables at the 
//      interface with linear interpolation
//
//      Parameters:
//          prim:   output parameter for the reconstructed primitive variables
// ============================================================================
void Interface::interpolatePrim(double * prim)
{
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

    if( fabs(this->normal.x) < 1.0e-12 )
        int i = 0;
}

// ============================================================================
//      This method computes the Values of the primary variables at the 
//      interface with pice wise constant reconstruction and averaging.
//
//      Parameters:
//          prim:   output parameter for the reconstructed primitive variables
// ============================================================================
__declspec(noinline) void Interface::reconstructPrimPiecewiseConstant(double * prim)
{
    // This method computes the Values of the primary variables at the interface
    // with pice wise constant reconstruction and averaging

    double distance = this->posDistance + this->negDistance;

    prim[0] = 0.5 * ( this->negCell->getPrim().rho 
                    + this->posCell->getPrim().rho );

    prim[1] = 0.5 * ( this->negCell->getPrim().U   
                    + this->posCell->getPrim().U   );

    prim[2] = 0.5 * ( this->negCell->getPrim().V   
                    + this->posCell->getPrim().V   );

    prim[3] = 0.5 * ( this->negCell->getPrim().L   
                    + this->posCell->getPrim().L   );

    if( fabs(this->normal.x) < 1.0e-12 )
        int i = 0;
}

// ============================================================================
//      This method computes the Values of the primary variables at the 
//      interface with picewise linear reconstruction based on the cell wise
//      gradients (computed by the least square approach) and averaging.
//
//      Parameters:
//          prim:   output parameter for the reconstructed primitive variables
// ============================================================================
void Interface::reconstructPrimPiecewiseLinear(double * prim)
{
    ConservedVariable negGradientX = this->negCell->getGradientX();
    ConservedVariable negGradientY = this->negCell->getGradientY();
    ConservedVariable posGradientX = this->posCell->getGradientX();
    ConservedVariable posGradientY = this->posCell->getGradientY();

    ConservedVariable cons;
    ConservedVariable negCons;
    ConservedVariable posCons;

    double negDx = this->center.x - this->negCell->getCenter().x;
    double negDy = this->center.y - this->negCell->getCenter().y;
    double posDx = this->center.x - this->posCell->getCenter().x;
    double posDy = this->center.y - this->posCell->getCenter().y;

    negCons.rho  = this->negCell->getCons().rho  + negGradientX.rho  * negDx + negGradientY.rho  * negDy;
    negCons.rhoU = this->negCell->getCons().rhoU + negGradientX.rhoU * negDx + negGradientY.rhoU * negDy;
    negCons.rhoV = this->negCell->getCons().rhoV + negGradientX.rhoV * negDx + negGradientY.rhoV * negDy;
    negCons.rhoE = this->negCell->getCons().rhoE + negGradientX.rhoE * negDx + negGradientY.rhoE * negDy;

    posCons.rho  = this->posCell->getCons().rho  + posGradientX.rho  * posDx + posGradientY.rho  * posDy;
    posCons.rhoU = this->posCell->getCons().rhoU + posGradientX.rhoU * posDx + posGradientY.rhoU * posDy;
    posCons.rhoV = this->posCell->getCons().rhoV + posGradientX.rhoV * posDx + posGradientY.rhoV * posDy;
    posCons.rhoE = this->posCell->getCons().rhoE + posGradientX.rhoE * posDx + posGradientY.rhoE * posDy;

    cons.rho  = 0.5 * ( negCons.rho  + posCons.rho  );
    cons.rhoU = 0.5 * ( negCons.rhoU + posCons.rhoU );
    cons.rhoV = 0.5 * ( negCons.rhoV + posCons.rhoV );
    cons.rhoE = 0.5 * ( negCons.rhoE + posCons.rhoE );

    prim[0] = this->cons2Prim(cons).rho;
    prim[1] = this->cons2Prim(cons).U;
    prim[2] = this->cons2Prim(cons).V;
    prim[3] = this->cons2Prim(cons).L;

    int i = 0;
}

// ====================================================================================================================
//
//                          Differentiation of conserved variables
//
// ====================================================================================================================

// ============================================================================
//      This method computes the normal derivatives of the conservative
//      Variables.
//      The derivatives are computed by two point central finite differences.
//
//      Parameters:
//          normalGradCons:   output parameters for the derivative
//          prim:             array of primitive variables
// ============================================================================
__declspec(noinline) void Interface::differentiateConsNormal(double* normalGradCons, double* prim)
{
    double distance = this->posDistance + this->negDistance;

    normalGradCons[0] = ( this->posCell->getCons().rho  - this->negCell->getCons().rho )  / ( distance * prim[0] );

    normalGradCons[1] = ( this->posCell->getCons().rhoU - this->negCell->getCons().rhoU ) / ( distance * prim[0] );

    normalGradCons[2] = ( this->posCell->getCons().rhoV - this->negCell->getCons().rhoV ) / ( distance * prim[0] );

    normalGradCons[3] = ( this->posCell->getCons().rhoE - this->negCell->getCons().rhoE ) / ( distance * prim[0] );
}

// ============================================================================
//      This method computes the normal derivatives of the conservative
//      Variables.
//      The derivatives are computed by three point central finite differences.
//      The values on the interface are assumed to be the average bet pos and 
//      neg side
//
//      Parameters:
//          normalGradCons:   output parameters for the derivative
//          prim:             array of primitive variables
// ============================================================================
void Interface::differentiateConsNormalThreePoint(double* normalGradCons, double * prim)
{
    double distanceFacctor = 0.5 * ( this->posDistance * this->posDistance + this->negDistance * this->negDistance ) 
                                 / ( this->posDistance * this->negDistance * ( this->posDistance + this->negDistance ) );

    normalGradCons[0] = distanceFacctor * ( this->posCell->getCons().rho  - this->negCell->getCons().rho )  / ( prim[0] );

    normalGradCons[1] = distanceFacctor * ( this->posCell->getCons().rhoU - this->negCell->getCons().rhoU ) / ( prim[0] );

    normalGradCons[2] = distanceFacctor * ( this->posCell->getCons().rhoV - this->negCell->getCons().rhoV ) / ( prim[0] );

    normalGradCons[3] = distanceFacctor * ( this->posCell->getCons().rhoE - this->negCell->getCons().rhoE ) / ( prim[0] );
}

// ============================================================================
//      This method computes the gradient of the conservative
//      Variables.
//      The derivatives are computed by averaging with decoupling correction.
//
//      Parameters:
//          normalGradCons:         output parameters for the normal derivative
//          tangentialGradCons:     output parameters for the tangential derivative
//          prim:                   array of primitive variables
// ============================================================================
void Interface::differentiateConsLeastSquare(double* normalGradCons, double* tangentialGradCons, double* prim)
{
    // ========================================================================
    // Interpolate Gradients
    // ========================================================================
    double distance = this->posDistance + this->negDistance;

    ConservedVariable interpolatedGradientX;
    ConservedVariable interpolatedGradientY;

    //interpolatedGradientX.rho  = ( this->posCell->getGradientX().rho  * this->negDistance + this->negCell->getGradientX().rho  * this->posDistance ) / ( distance );
    //interpolatedGradientX.rhoU = ( this->posCell->getGradientX().rhoU * this->negDistance + this->negCell->getGradientX().rhoU * this->posDistance ) / ( distance );
    //interpolatedGradientX.rhoV = ( this->posCell->getGradientX().rhoV * this->negDistance + this->negCell->getGradientX().rhoV * this->posDistance ) / ( distance );
    //interpolatedGradientX.rhoE = ( this->posCell->getGradientX().rhoE * this->negDistance + this->negCell->getGradientX().rhoE * this->posDistance ) / ( distance );

    //interpolatedGradientY.rho  = ( this->posCell->getGradientY().rho  * this->negDistance + this->negCell->getGradientY().rho  * this->posDistance ) / ( distance );
    //interpolatedGradientY.rhoU = ( this->posCell->getGradientY().rhoU * this->negDistance + this->negCell->getGradientY().rhoU * this->posDistance ) / ( distance );
    //interpolatedGradientY.rhoV = ( this->posCell->getGradientY().rhoV * this->negDistance + this->negCell->getGradientY().rhoV * this->posDistance ) / ( distance );
    //interpolatedGradientY.rhoE = ( this->posCell->getGradientY().rhoE * this->negDistance + this->negCell->getGradientY().rhoE * this->posDistance ) / ( distance );

    interpolatedGradientX.rho  = 0.5 * ( this->posCell->getGradientX().rho  + this->negCell->getGradientX().rho  );
    interpolatedGradientX.rhoU = 0.5 * ( this->posCell->getGradientX().rhoU + this->negCell->getGradientX().rhoU );
    interpolatedGradientX.rhoV = 0.5 * ( this->posCell->getGradientX().rhoV + this->negCell->getGradientX().rhoV );
    interpolatedGradientX.rhoE = 0.5 * ( this->posCell->getGradientX().rhoE + this->negCell->getGradientX().rhoE );
    
    interpolatedGradientY.rho  = 0.5 * ( this->posCell->getGradientY().rho  + this->negCell->getGradientY().rho  );
    interpolatedGradientY.rhoU = 0.5 * ( this->posCell->getGradientY().rhoU + this->negCell->getGradientY().rhoU );
    interpolatedGradientY.rhoV = 0.5 * ( this->posCell->getGradientY().rhoV + this->negCell->getGradientY().rhoV );
    interpolatedGradientY.rhoE = 0.5 * ( this->posCell->getGradientY().rhoE + this->negCell->getGradientY().rhoE );
    // ========================================================================
    
    // ========================================================================
    // Decoupling correction as given in Blazeks Book
    // ========================================================================
    double dx = this->posCell->getCenter().x - this->negCell->getCenter().x;
    double dy = this->posCell->getCenter().y - this->negCell->getCenter().y;

    // Eq. 5.71 in Blazek
    ConservedVariable directionalGradient;
    directionalGradient.rho  = ( this->posCell->getCons().rho  - this->negCell->getCons().rho  ) / sqrt( dx*dx + dy*dy );
    directionalGradient.rhoU = ( this->posCell->getCons().rhoU - this->negCell->getCons().rhoU ) / sqrt( dx*dx + dy*dy );
    directionalGradient.rhoV = ( this->posCell->getCons().rhoV - this->negCell->getCons().rhoV ) / sqrt( dx*dx + dy*dy );
    directionalGradient.rhoE = ( this->posCell->getCons().rhoE - this->negCell->getCons().rhoE ) / sqrt( dx*dx + dy*dy );

    // Eq. 5.72 in Blazek
    dx /= sqrt( dx*dx + dy*dy );
    dy /= sqrt( dx*dx + dy*dy );

    // eq. 5.73 in Blazek
    interpolatedGradientX.rho  -= ( interpolatedGradientX.rho  * dx + interpolatedGradientY.rho  * dy - directionalGradient.rho  ) * dx;
    interpolatedGradientX.rhoU -= ( interpolatedGradientX.rhoU * dx + interpolatedGradientY.rhoU * dy - directionalGradient.rhoU ) * dx;
    interpolatedGradientX.rhoV -= ( interpolatedGradientX.rhoV * dx + interpolatedGradientY.rhoV * dy - directionalGradient.rhoV ) * dx;
    interpolatedGradientX.rhoE -= ( interpolatedGradientX.rhoE * dx + interpolatedGradientY.rhoE * dy - directionalGradient.rhoE ) * dx;

    interpolatedGradientY.rho  -= ( interpolatedGradientX.rho  * dx + interpolatedGradientY.rho  * dy - directionalGradient.rho  ) * dy;
    interpolatedGradientY.rhoU -= ( interpolatedGradientX.rhoU * dx + interpolatedGradientY.rhoU * dy - directionalGradient.rhoU ) * dy;
    interpolatedGradientY.rhoV -= ( interpolatedGradientX.rhoV * dx + interpolatedGradientY.rhoV * dy - directionalGradient.rhoV ) * dy;
    interpolatedGradientY.rhoE -= ( interpolatedGradientX.rhoE * dx + interpolatedGradientY.rhoE * dy - directionalGradient.rhoE ) * dy;
    // ========================================================================
    
    // ========================================================================
    // transformation from global into local coordinatesystem and normalization
    //    by projection onto normal and tangential vectors
    // ========================================================================
    normalGradCons[0]     = (   this->normal.x * interpolatedGradientX.rho  + this->normal.y * interpolatedGradientY.rho  ) / prim[0];
    normalGradCons[1]     = (   this->normal.x * interpolatedGradientX.rhoU + this->normal.y * interpolatedGradientY.rhoU ) / prim[0];
    normalGradCons[2]     = (   this->normal.x * interpolatedGradientX.rhoV + this->normal.y * interpolatedGradientY.rhoV ) / prim[0];
    normalGradCons[3]     = (   this->normal.x * interpolatedGradientX.rhoE + this->normal.y * interpolatedGradientY.rhoE ) / prim[0];

    tangentialGradCons[0] = ( - this->normal.y * interpolatedGradientX.rho  + this->normal.x * interpolatedGradientY.rho  ) / prim[0];
    tangentialGradCons[1] = ( - this->normal.y * interpolatedGradientX.rhoU + this->normal.x * interpolatedGradientY.rhoU ) / prim[0];
    tangentialGradCons[2] = ( - this->normal.y * interpolatedGradientX.rhoV + this->normal.x * interpolatedGradientY.rhoV ) / prim[0];
    tangentialGradCons[3] = ( - this->normal.y * interpolatedGradientX.rhoE + this->normal.x * interpolatedGradientY.rhoE ) / prim[0];
    // ========================================================================

    int breakPoint = 0;
}

// ====================================================================================================================
//
//                          Flux computation
//
// ====================================================================================================================

// ============================================================================
//      This method computes the Fluxes over an Interface
//
//      Parameters:
//          dt:   time step
// ============================================================================
void Interface::computeInternalFlux(double dt)
{
    const int NUMBER_OF_MOMENTS = 7;

    double prim[4];
    double primTest[4] = {0.0, 0.0, 0.0, 0.0};

    double normalGradCons[4];
    double normalGradConsTest[4];

    double tangentialGradCons[4] = {0.0, 0.0, 0.0, 0.0};

    double timeGrad[4];

    double a[4];
    double b[4] = {0.0, 0.0, 0.0, 0.0};
    double A[4];

    double MomentU[NUMBER_OF_MOMENTS];
    double MomentV[NUMBER_OF_MOMENTS];
    double MomentXi[NUMBER_OF_MOMENTS];

    // ========================================================================
    //          interpolated primitive variables at the interface
    // ========================================================================
    //this->interpolatePrim(prim);
    this->reconstructPrimPiecewiseConstant(prim);
    //this->reconstructPrimPiecewiseConstant(primTest);
    //this->reconstructPrimPiecewiseLinear(prim);
    // ========================================================================
    
    // ========================================================================
    //          compute spacial gradients of the conservative varibles
    // ========================================================================
    this->differentiateConsNormal(normalGradCons, prim);
    //this->differentiateConsNormal(normalGradConsTest, prim);
    //this->differentiateConsNormalThreePoint(normalGradCons, prim);
    //this->differentiateConsLeastSquare(normalGradCons, tangentialGradCons, prim);
    // ========================================================================
    
    // ========================================================================
    //          Some Tests
    // ========================================================================
    //PrimitiveVariable posPrim = this->posCell->getPrim();
    //PrimitiveVariable negPrim = this->negCell->getPrim();

    ////if( fabs(this->normal.x - 1.0) > 1.0e-12 )
    //if(this->ID == 9)
    //    int breakPoint = 0;

    //this->transformGlobal2Local( (double*)&posPrim );
    //this->transformGlobal2Local( (double*)&negPrim );
    // ========================================================================
    
    // ========================================================================
    //          Momentum Transformation in local coordinate system
    // ========================================================================
    transformGlobal2Local(prim);
    transformGlobal2Local(normalGradCons);
    transformGlobal2Local(tangentialGradCons);
    // ========================================================================
    
    // ========================================================================
    //          compute spacial micro slopes
    //              a = a1 + a2 u + a3 v + 0.5 a4 (u^2 + v^2 + xi^2)
    // ========================================================================
    this->computeMicroSlope(prim, normalGradCons, a);
    //this->computeMicroSlope(prim, tangentialGradCons, b);
    // ========================================================================
    
    // ========================================================================
    //          comoute moments of the equilibrium distribution
    // ========================================================================
    this->computeMoments(prim, MomentU, MomentV, MomentXi, NUMBER_OF_MOMENTS);
    // ========================================================================

    // ========================================================================
    //          compute time derivative and temporal micro slopes
    //              A = A1 + A2 u + A3 v + 0.5 A4 (u^2 + v^2 + xi^2)
    // ========================================================================
    this->computeTimeDerivative(prim, MomentU, MomentV, MomentXi, a, b, timeGrad);

    this->computeMicroSlope(prim, timeGrad, A);
    // ========================================================================

    // ========================================================================
    // Relaxation time as in the Rayleigh-Bernard-Paper (Xu, Lui, 1999)
    // ========================================================================
    double tau = 2.0*prim[3] * this->fluidParam.nu;
    // ========================================================================
    
    // ========================================================================
    //          compute time integration Coefficients
    // ========================================================================
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };
    // ========================================================================

    // ========================================================================
    //          compute mass and momentum fluxes
    // ========================================================================
    this->assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients, prim, tau);
    // ========================================================================

    if(this->ID == 9)
        int breakPoint = 1;
    
    // ========================================================================
    //          transform momentum Flux components back to global system
    // ========================================================================
    transformLocal2Global(this->timeIntegratedFlux);
    transformLocal2Global(this->timeIntegratedFlux_1);
    transformLocal2Global(this->timeIntegratedFlux_2);
    transformLocal2Global(this->timeIntegratedFlux_3);
    transformLocal2Global(this->FluxDensity);
    // ========================================================================
    
    if(this->ID == 9)
        int breakPoint = 1;
}

// ============================================================================
//      This method computes the moments of the equilibrium distribution.
//      The formulas can be found in (Xu, 2001).
//
//      Parameters:
//          prim:           input parameter for primitive variables
//          MomentU:        output parameter for the moments with respect to u
//          MomentV:        output parameter for the moments with respect to v
//          MomentXi:       output parameter for the moments with respect to xi
//          numberMoments:  input parameter for how many moments shall be computed
// ============================================================================
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

// ====================================================================================================================
//
//                          Data Manipulation
//
// ====================================================================================================================

// ============================================================================
//      This method transforms the second and third element in vec (vec[1] and
//      vec[2]) into the local frame of reference.
//      This is done by projection onto normal and tangential vectors.
//
//      Parameters:
//          vec:   array in which the second and third element are transformed
// ============================================================================
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

// ============================================================================
//      This method transforms the second and third element in vec (vec[1] and
//      vec[2]) into the global frame of reference.
//      This is done by the tanspoed projection onto normal and tangential
//      vectors.
//
//      Parameters:
//          vec:   array in which the second and third element are transformed
// ============================================================================
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

// ============================================================================
//      This method returns a set of primitive variables equivalent to the
//      conservative variables in cons.
//
//      Parameters:
//          cons:   conserved variable
// ============================================================================
PrimitiveVariable Interface::cons2Prim(ConservedVariable cons)
{
    PrimitiveVariable prim;
    prim.rho = cons.rho;
    prim.U   = cons.rhoU / cons.rho;
    prim.V   = cons.rhoV / cons.rho;

    prim.L = ( this->fluidParam.K + 2.0 ) * cons.rho / ( 4.0 * ( cons.rhoE - 0.5 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV ) / cons.rho ) );

    return prim;
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          get Methods
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ====================================================================================================================
//
//                          get Connectivity
//
// ====================================================================================================================

idType Interface::getID()
{
    return this->ID;
}

// ============================================================================
//      This method returns the i-th node of the interface
// ============================================================================
Node * Interface::getNode(int i)
{
    return this->nodes[i];
}

Cell * Interface::getCell(idType i)
{
    if(i == 0) return this->posCell;
    else       return this->negCell;
}

// ============================================================================
//      This method returns the opposite cell of the asking cell at this 
//      interface.
// ============================================================================
Cell * Interface::getNeigborCell(Cell * askingCell)
{
    if (posCell == askingCell)
        return negCell;
    else
        return posCell;
}

// ============================================================================
//      This method is intended for boundary interfaces to get the cell that 
//      is in the domain.
// ============================================================================
Cell * Interface::getCellInDomain()
{
    if( !this->isBoundaryInterface() )
        return NULL;
    if ( !posCell->isGhostCell() )
        return posCell;
    if ( !negCell->isGhostCell() )
        return negCell;
}

// ============================================================================
//      This method is intended for periodic boundary interfaces. It gets the
//      corresponding cell the opposite side of the domain.
// ============================================================================
Cell * Interface::getPeriodicCell()
{
    Cell* tmp = NULL;

    if( fabs( posDistance - this->distance( this->posCell->getCenter() ) ) > 1.0e-12  )
    {
        tmp = posCell;
        posCell = NULL;
    }

    if( fabs( negDistance - this->distance( this->negCell->getCenter() ) ) > 1.0e-12  )
    {
        tmp = negCell;
        negCell = NULL;
    }

    return tmp;
}

// ============================================================================
//      This method returns a vector from the interface center to the center
//      of the adjacent cell on the positive side.
// ============================================================================
Node Interface::getPosConnectivity()
{
    Node vector;
    vector.x = this->posCell->getCenter().x - this->center.x;
    vector.y = this->posCell->getCenter().y - this->center.y;
    return vector;
}

// ============================================================================
//      This method returns a vector from the interface center to the center
//      of the adjacent cell on the negative side.
// ============================================================================
Node Interface::getNegConnectivity()
{
    Node vector;
    vector.x = this->negCell->getCenter().x - this->center.x;
    vector.y = this->negCell->getCenter().y - this->center.y;
    return vector;
}

// ============================================================================
//      This method returns a pointer to the boundary condition of this 
//      interface.
// ============================================================================
BoundaryCondition * Interface::getBoundaryCondition()
{
    return this->BoundaryConditionPointer;
}

// ============================================================================
//      This method returns, wheter this interface connects two ghostcells.
// ============================================================================
bool Interface::isGhostInterface()
{
    if (this->isBoundaryInterface()) return false;

    return this->posCell->isGhostCell() && this->negCell->isGhostCell();
}

// ============================================================================
//      This method returns, wheter this interface is a boundary itnerface.
//      It does not recognize periodic boundaries without ghost cells.
// ============================================================================
bool Interface::isBoundaryInterface()
{
    return this->BoundaryConditionPointer != NULL;
}

// ====================================================================================================================
//
//                          get Geometry
//
// ====================================================================================================================

// ============================================================================
//      This method returns the normal vector of the interface.
// ============================================================================
Node Interface::getNormal()
{
    return normal;
}

// ============================================================================
//      This method returns the center of this interface.
// ============================================================================
Node Interface::getCenter()
{
    return this->center;
}

// ============================================================================
//      This method is inteded for the construction of ghost cells. it returns
//      a normal vector pointing out of the domain with a length that is equal
//      to the distance between the interface center and the adjacent cell
//      Center in the domain.
// ============================================================================
Node Interface::getScaledNormal()
{
    if(this->posCell == NULL)
        return Node(   this->normal.x * this->distance(negCell->getCenter()),   this->normal.y * this->distance(negCell->getCenter()) );
    else
        return Node( - this->normal.x * this->distance(posCell->getCenter()), - this->normal.y * this->distance(posCell->getCenter()) );
}

// ============================================================================
//      This method returns the interface area.
// ============================================================================
double Interface::getArea()
{
    return this->area;
}

double Interface::getDistance2CellCenter(idType i)
{
    if(i == 0) return this->posDistance;
    else       return this->negDistance;
}

// ============================================================================
//      This method projects the vector from the interface center to some point
//      and returns its length.
// ============================================================================
double Interface::distance(Node point)
{
    // Compute the projected distance of a point on the normal of the interface
    return fabs( ( this->center.x - point.x )*( this->normal.x )
               + ( this->center.y - point.y )*( this->normal.y ) );
}

// ====================================================================================================================
//
//                          get Data
//
// ====================================================================================================================

// ============================================================================
//      This method returns the time integrated Flux over this interface.
// ============================================================================
ConservedVariable Interface::getTimeIntegratedFlux()
{
    ConservedVariable tmp;
    tmp.rho  = this->timeIntegratedFlux[0];
    tmp.rhoU = this->timeIntegratedFlux[1];
    tmp.rhoV = this->timeIntegratedFlux[2];
    tmp.rhoE = this->timeIntegratedFlux[3];
    return tmp;
}

// ============================================================================
//      This method returns the time integrated Flux over this interface,
//      that is related to convetive transport.
// ============================================================================
ConservedVariable Interface::getTimeIntegratedFlux_1()
{
    ConservedVariable tmp;
    tmp.rho  = this->timeIntegratedFlux_1[0];
    tmp.rhoU = this->timeIntegratedFlux_1[1];
    tmp.rhoV = this->timeIntegratedFlux_1[2];
    tmp.rhoE = this->timeIntegratedFlux_1[3];
    return tmp;
}

// ============================================================================
//      This method returns the time integrated Flux over this interface,
//      that is related transport due to spacial gradients.
// ============================================================================
ConservedVariable Interface::getTimeIntegratedFlux_2()
{
    ConservedVariable tmp;
    tmp.rho  = this->timeIntegratedFlux_2[0];
    tmp.rhoU = this->timeIntegratedFlux_2[1];
    tmp.rhoV = this->timeIntegratedFlux_2[2];
    tmp.rhoE = this->timeIntegratedFlux_2[3];
    return tmp;
}

// ============================================================================
//      This method returns the time integrated Flux over this interface,
//      that is related transport due to temporal gradients.
// ============================================================================
ConservedVariable Interface::getTimeIntegratedFlux_3()
{
    ConservedVariable tmp;
    tmp.rho  = this->timeIntegratedFlux_3[0];
    tmp.rhoU = this->timeIntegratedFlux_3[1];
    tmp.rhoV = this->timeIntegratedFlux_3[2];
    tmp.rhoE = this->timeIntegratedFlux_3[3];
    return tmp;
}

// ============================================================================
//      This method returns the Flux density over this interface.
// ============================================================================
ConservedVariable Interface::getFluxDensity()
{
    ConservedVariable tmp;
    tmp.rho =  this->FluxDensity[0];
    tmp.rhoU = this->FluxDensity[1];
    tmp.rhoV = this->FluxDensity[2];
    tmp.rhoE = this->FluxDensity[3];
    return tmp;
}

// ============================================================================
//      This method returns 1.0 or -1.0 depending on whether the asking cell
//      is on the positive or negative side of the interface
// ============================================================================
double Interface::getFluxSign(Cell * askingCell)
{
    double sign = ( askingCell->getCenter().x - this->center.x ) * this->normal.x
                + ( askingCell->getCenter().y - this->center.y ) * this->normal.y ;

    return sign/fabs(sign);
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          output Methods
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method returns a string containing some information on this 
//      interface.
// ============================================================================
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

// ============================================================================
//      This method writes the interface Center to a string.
// ============================================================================
string Interface::writeCenter()
{
    ostringstream tmp;
    tmp << this->center.x << " " << this->center.y << " 0.0\n";
    return tmp.str();
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
