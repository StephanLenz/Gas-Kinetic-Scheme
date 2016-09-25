
#include "Cell.h"
#include "Interface.h"
#include <sstream>
#include <iostream>
#include <algorithm>    // min()

using namespace std;

Cell::Cell()
{
    memset(InterfaceList, NULL, 4 * sizeof(Interface*));
}

Cell::Cell(InterfaceType interfaceType, float2** nodes, BoundaryCondition* BC, FluidParameter fluidParam)
{
    // ========================================================================
    //                  Copy attributes
    // ========================================================================
    memset(InterfaceList, NULL, 4 * sizeof(Interface*));

    this->nInterfaces = 0;

    for(int i = 0; i < 4; i++)
        this->nodes[i] = nodes[i];

    this->BoundaryContitionPointer = BC;

    this->fluidParam = fluidParam;

    this->interfaceType = interfaceType;
    // ========================================================================
    
    // ========================================================================
    //                  Compute Centroid and Volume of the Quad
    // ========================================================================
    float2 triCenter[2];
    triCenter[0].x =  (this->nodes[0]->x + this->nodes[1]->x +                     this->nodes[3]->x) / 3.0;
    triCenter[0].y =  (this->nodes[0]->y + this->nodes[1]->y +                     this->nodes[3]->y) / 3.0;
    triCenter[1].x =  (                    this->nodes[1]->x + this->nodes[2]->x + this->nodes[3]->x) / 3.0;
    triCenter[1].y =  (                    this->nodes[1]->y + this->nodes[2]->y + this->nodes[3]->y) / 3.0;

    double triVolume[2];
    triVolume[0] = 0.5 * fabs( this->nodes[0]->x * ( this->nodes[1]->y - this->nodes[3]->y ) 
                             + this->nodes[1]->x * ( this->nodes[3]->y - this->nodes[0]->y ) 
                             + this->nodes[3]->x * ( this->nodes[0]->y - this->nodes[1]->y ) );
    triVolume[1] = 0.5 * fabs( this->nodes[2]->x * ( this->nodes[3]->y - this->nodes[1]->y ) 
                             + this->nodes[3]->x * ( this->nodes[1]->y - this->nodes[2]->y ) 
                             + this->nodes[1]->x * ( this->nodes[2]->y - this->nodes[3]->y ) );

    this->volume = triVolume[0] + triVolume[1];
    this->center.x = ( triCenter[0].x * triVolume[0] + triCenter[1].x * triVolume[1] ) / this->volume;
    this->center.y = ( triCenter[0].y * triVolume[0] + triCenter[1].y * triVolume[1] ) / this->volume;
    // ========================================================================
    
    this->updateVal.rho  = 0.0;
    this->updateVal.rhoU = 0.0;
    this->updateVal.rhoV = 0.0;
    this->updateVal.rhoE = 0.0;

    for(int i = 0; i < 4; i++)
    {
        ((double*)&this->gradientX)[i] = 0.0;
        ((double*)&this->gradientY)[i] = 0.0;
    }

    // minDx must be computed in annother step, when the interfaces are created and connected
    this->minDx = 1.0e99;
}


Cell::~Cell()
{
}

void Cell::update(double dt)
{
    // ========================================================================
    //                      Update conservative Variables
    // ========================================================================
    this->cons[0] += ( this->InterfaceList[0]->getFluxSign(this) * this->InterfaceList[0]->getTimeIntegratedFlux().rho
                     + this->InterfaceList[1]->getFluxSign(this) * this->InterfaceList[1]->getTimeIntegratedFlux().rho
                     + this->InterfaceList[2]->getFluxSign(this) * this->InterfaceList[2]->getTimeIntegratedFlux().rho
                     + this->InterfaceList[3]->getFluxSign(this) * this->InterfaceList[3]->getTimeIntegratedFlux().rho
                     ) / this->volume;

    this->cons[1] += ( this->InterfaceList[0]->getFluxSign(this) * this->InterfaceList[0]->getTimeIntegratedFlux().rhoU
                     + this->InterfaceList[1]->getFluxSign(this) * this->InterfaceList[1]->getTimeIntegratedFlux().rhoU
                     + this->InterfaceList[2]->getFluxSign(this) * this->InterfaceList[2]->getTimeIntegratedFlux().rhoU
                     + this->InterfaceList[3]->getFluxSign(this) * this->InterfaceList[3]->getTimeIntegratedFlux().rhoU
                     ) / this->volume;

    this->cons[2] += ( this->InterfaceList[0]->getFluxSign(this) * this->InterfaceList[0]->getTimeIntegratedFlux().rhoV
                     + this->InterfaceList[1]->getFluxSign(this) * this->InterfaceList[1]->getTimeIntegratedFlux().rhoV
                     + this->InterfaceList[2]->getFluxSign(this) * this->InterfaceList[2]->getTimeIntegratedFlux().rhoV
                     + this->InterfaceList[3]->getFluxSign(this) * this->InterfaceList[3]->getTimeIntegratedFlux().rhoV
                     ) / this->volume;

    this->cons[3] += ( this->InterfaceList[0]->getFluxSign(this) * this->InterfaceList[0]->getTimeIntegratedFlux().rhoE
                     + this->InterfaceList[1]->getFluxSign(this) * this->InterfaceList[1]->getTimeIntegratedFlux().rhoE
                     + this->InterfaceList[2]->getFluxSign(this) * this->InterfaceList[2]->getTimeIntegratedFlux().rhoE
                     + this->InterfaceList[3]->getFluxSign(this) * this->InterfaceList[3]->getTimeIntegratedFlux().rhoE
                     ) / this->volume;
    // ========================================================================

    // ========================================================================
    //                        Some Test
    // ========================================================================
    ConservedVariable Flux0 = this->InterfaceList[0]->getTimeIntegratedFlux();
    ConservedVariable Flux1 = this->InterfaceList[1]->getTimeIntegratedFlux();
    ConservedVariable Flux2 = this->InterfaceList[2]->getTimeIntegratedFlux();
    ConservedVariable Flux3 = this->InterfaceList[3]->getTimeIntegratedFlux();

    double Sign0 = this->InterfaceList[0]->getFluxSign(this); 
    double Sign1 = this->InterfaceList[1]->getFluxSign(this); 
    double Sign2 = this->InterfaceList[2]->getFluxSign(this); 
    double Sign3 = this->InterfaceList[3]->getFluxSign(this); 

    ConservedVariable update;

    update.rho     = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux().rho
                     + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux().rho
                     + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux().rho
                     + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux().rho
                     ) / this->volume;

    update.rhoU    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux().rhoU
                     + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux().rhoU
                     + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux().rhoU
                     + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux().rhoU
                     ) / this->volume;

    update.rhoV    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux().rhoV
                     + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux().rhoV
                     + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux().rhoV
                     + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux().rhoV
                     ) / this->volume;

    update.rhoE    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux().rhoE
                     + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux().rhoE
                     + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux().rhoE
                     + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux().rhoE
                     ) / this->volume;

    //// ========================================================================

    //ConservedVariable update_1;

    //update_1.rho     = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_1().rho
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_1().rho
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_1().rho
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_1().rho
    //                 ) / this->volume;

    //update_1.rhoU    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_1().rhoU
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_1().rhoU
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_1().rhoU
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_1().rhoU
    //                 ) / this->volume;

    //update_1.rhoV    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_1().rhoV
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_1().rhoV
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_1().rhoV
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_1().rhoV
    //                 ) / this->volume;

    //update_1.rhoE    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_1().rhoE
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_1().rhoE
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_1().rhoE
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_1().rhoE
    //                 ) / this->volume;

    //ConservedVariable update_2;

    //update_2.rho     = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_2().rho
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_2().rho
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_2().rho
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_2().rho
    //                 ) / this->volume;

    //update_2.rhoU    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_2().rhoU
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_2().rhoU
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_2().rhoU
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_2().rhoU
    //                 ) / this->volume;

    //update_2.rhoV    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_2().rhoV
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_2().rhoV
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_2().rhoV
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_2().rhoV
    //                 ) / this->volume;

    //update_2.rhoE    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_2().rhoE
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_2().rhoE
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_2().rhoE
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_2().rhoE
    //                 ) / this->volume;

    //ConservedVariable update_3;

    //update_3.rho     = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_3().rho
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_3().rho
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_3().rho
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_3().rho
    //                 ) / this->volume;

    //update_3.rhoU    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_3().rhoU
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_3().rhoU
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_3().rhoU
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_3().rhoU
    //                 ) / this->volume;

    //update_3.rhoV    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_3().rhoV
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_3().rhoV
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_3().rhoV
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_3().rhoV
    //                 ) / this->volume;

    //update_3.rhoE    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_3().rhoE
    //                   + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_3().rhoE
    //                   + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_3().rhoE
    //                   + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_3().rhoE
    //                 ) / this->volume;

    //ConservedVariable update_test;

    //update_test.rho  = update_1.rho  + update_2.rho  + update_3.rho;
    //update_test.rhoU = update_1.rhoU + update_2.rhoU + update_3.rhoU;
    //update_test.rhoV = update_1.rhoV + update_2.rhoV + update_3.rhoV;
    //update_test.rhoE = update_1.rhoE + update_2.rhoE + update_3.rhoE;

    //if(update_test.rho  < 1.0e-10 && fabs(this->cons[0] - this->cons_old[0]) < 1.0e-10) update_test.rho  = 0.0;
    //if(update_test.rhoU < 1.0e-10 && fabs(this->cons[1] - this->cons_old[1]) < 1.0e-10) update_test.rhoU = 0.0;
    //if(update_test.rhoV < 1.0e-10 && fabs(this->cons[2] - this->cons_old[2]) < 1.0e-10) update_test.rhoV = 0.0;
    //if(update_test.rhoE < 1.0e-10 && fabs(this->cons[3] - this->cons_old[3]) < 1.0e-10) update_test.rhoE = 0.0;

    // ========================================================================

    //this->cons[0] += this->updateVal.rho  / this->volume;
    //this->cons[1] += this->updateVal.rhoU / this->volume;
    //this->cons[2] += this->updateVal.rhoV / this->volume;
    //this->cons[3] += this->updateVal.rhoE / this->volume;

    //int l = 0;

    this->updateVal.rho  = update.rho ;
    this->updateVal.rhoU = update.rhoU;
    this->updateVal.rhoV = update.rhoV;
    this->updateVal.rhoE = update.rhoE;

    int i = 0;
    // ========================================================================

    // ========================================================================
    //                        Compute Residual and store new values
    // ========================================================================
    this->residual.rho  = fabs(this->cons[0] - this->cons_old[0]);
    this->residual.rhoU = fabs(this->cons[1] - this->cons_old[1]);
    this->residual.rhoV = fabs(this->cons[2] - this->cons_old[2]);
    this->residual.rhoE = fabs(this->cons[3] - this->cons_old[3]);

    int k = 0;

    // store values of this time step for residual computation in the next one
    this->cons_old[0] = this->cons[0];
    this->cons_old[1] = this->cons[1];
    this->cons_old[2] = this->cons[2];
    this->cons_old[3] = this->cons[3];
    // ========================================================================

    int j = 0;

}

void Cell::applyBoundaryCondition()
{

    PrimitiveVariable primNeighbor = this->findNeighborInDomain()->getPrim();
    PrimitiveVariable prim;

    switch ( this->BoundaryContitionPointer->getType() )
    {
        // ====================================================================
        case wall:
        {
            PrimitiveVariable boundaryValue = this->BoundaryContitionPointer->getValue();

            prim.rho = primNeighbor.rho;
            prim.U   = 2.0*boundaryValue.U - primNeighbor.U;
            prim.V   = 2.0*boundaryValue.V - primNeighbor.V;
            prim.L   = primNeighbor.L;

            break;
        }
        // ====================================================================
        case isothermalWall:
        {
            PrimitiveVariable boundaryValue = this->BoundaryContitionPointer->getValue();

            prim.rho = primNeighbor.rho * ( 2.0*boundaryValue.L - primNeighbor.L ) / primNeighbor.L;
            prim.U   = 2.0*boundaryValue.U - primNeighbor.U;
            prim.V   = 2.0*boundaryValue.V - primNeighbor.V;
            prim.L   = 2.0*boundaryValue.L - primNeighbor.L;
                
            break;
        }
    }
    
    this->computeCons(prim);
}

void Cell::applyForcing(double dt)
{
    float2 Force;
    
    // ========================================================================
    //                  Compute Primitive Variables (esp. Temp before forcing)
    // ========================================================================
    PrimitiveVariable prim = this->getPrim();

    // ========================================================================
    //                  Compute Volumeforce from acceleration
    // ========================================================================
    Force.x = this->fluidParam.Force.x              *   this->cons[0]
            + this->fluidParam.BoussinesqForce.x    * ( this->cons[0] - this->fluidParam.rhoReference );
    Force.y = this->fluidParam.Force.y              *   this->cons[0]
            + this->fluidParam.BoussinesqForce.y    * ( this->cons[0] - this->fluidParam.rhoReference );
    // ========================================================================
    
    // ========================================================================
    //                  Apply Forcing to momentum components
    // ========================================================================
    this->cons[1] += dt * Force.x;
    this->cons[2] += dt * Force.y;
    // ========================================================================

    // ========== some Tests ==================================================
    double xForce = dt * Force.x;
    double yForce = dt * Force.y;

    int i = 0;
    // ========================================================================

    // ========================================================================
    //                  compute new Energy with increased momentum
    // ========================================================================
    this->cons[3] = prim.rho * (this->fluidParam.K + 2.0) / (4.0*prim.L)
                  + 0.5 * ( this->cons[1] * this->cons[1] 
                          + this->cons[2] * this->cons[2]
                          ) / prim.rho;
    // ========================================================================
    
}

void Cell::addInterface(Interface* newInterface)
{
	this->InterfaceList[this->nInterfaces++] = newInterface;
}

void Cell::computeMinDx()
{
    // ========================================================================
    //                  Compute minimal length of the Quad
    // ========================================================================
    if(! this->isGhostCell() )
    {
        for(int i = 0; i < 4; i++)          // loop over interfaces
        {
            for( int j = 0; j < 4; j++)     // loop over nodes not on the interface
            {
                if(  this->nodes[j] != this->InterfaceList[i]->getNode(0) 
                  && this->nodes[j] != this->InterfaceList[i]->getNode(1) )
                {
                    float2 node( this->nodes[j]->x, this->nodes[j]->y );
                    this->minDx = min( this->minDx, this->InterfaceList[i]->distance(node) );
                }
            }
        }
    }
    int i = 0;
}

void Cell::setValues(double rho, double u, double v, double L)
{
    PrimitiveVariable prim;
	prim.rho = rho;
	prim.U   = u;
	prim.V   = v;
	prim.L   = L;

    this->computeCons(prim);

    this->cons_old[0] = this->cons[0];
    this->cons_old[1] = this->cons[1];
    this->cons_old[2] = this->cons[2];
    this->cons_old[3] = this->cons[3];
}

void Cell::computeCons(PrimitiveVariable prim)
{
    this->cons[0] = prim.rho;
    this->cons[1] = prim.rho * prim.U; 
    this->cons[2] = prim.rho * prim.V;
    // inverse of eq. in GKS Book page 79 at the bottom
    this->cons[3] = prim.rho * (this->fluidParam.K + 2.0) / (4.0*prim.L)
                  + 0.5 * ( this->cons[1] * this->cons[1] 
                          + this->cons[2] * this->cons[2]
                          )/prim.rho;
}

double Cell::getLocalTimestep()
{
    PrimitiveVariable prim = this->getPrim();

    double U_max = sqrt(prim.U*prim.U + prim.V*prim.V);
    double c_s   = sqrt( 1.0 / ( 2.0*prim.L ) );                      // c_s = sqrt(RT) = c_s = sqrt(1/2lambda)

    double localTimestep = minDx / ( U_max + c_s + 2.0*this->fluidParam.nu / minDx );

    // ========================================================================
    
    return localTimestep;
}

void Cell::addFlux(double * Flux, double sign, Interface* origin)
{
    int  j = 0;
    bool flag = false;

    for(int i = 0; i < 4; ++i)
    {
        flag = flag || this->InterfaceList[i] == origin;
    }
    
    if(flag)
    {
        this->updateVal.rho  += sign * Flux[0];
        this->updateVal.rhoU += sign * Flux[1];
        this->updateVal.rhoV += sign * Flux[2];
        this->updateVal.rhoE += sign * Flux[3];

        int i = 0;
    }

    int i = 0;
}

void Cell::computeLeastSquareCoefficients()
{
    this->r11 = 0.0;
    this->r12 = 0.0;
    this->r22 = 0.0;

    for(int i = 0; i < 4; i++)
    {
        double dx = this->InterfaceList[i]->getNeigborCell(this)->getCenter().x - this->center.x;
        double dy = this->InterfaceList[i]->getNeigborCell(this)->getCenter().y - this->center.y;

        this->r11 += dx*dx;
        this->r12 += dx*dy;
        this->r22 += dy*dy;
    }

    this->r11 = sqrt(this->r11);
    this->r12 = this->r12 / r11;
    this->r22 = sqrt(this->r22 - r12*r12);
}

void Cell::computeGradients()
{
    // loop over conserved variables
    for(int i = 0; i < 4; i++)
    {
        ((double*)&this->gradientX)[i] = 0.0;
        ((double*)&this->gradientY)[i] = 0.0;

        for(int j = 0; j < 4; j++)
        {
            double dx = this->InterfaceList[j]->getNeigborCell(this)->getCenter().x - this->center.x;
            double dy = this->InterfaceList[j]->getNeigborCell(this)->getCenter().y - this->center.y;

            double dW = ((double*)&this->InterfaceList[j]->getNeigborCell(this)->getCons())[i] - this->cons[i];

            double w1 = dx / (r11*r11) - r12/r11 * ( dy - r12/r11 * dx ) / (r22*r22);
            double w2 =                            ( dy - r12/r11 * dx ) / (r22*r22);

            ((double*)&this->gradientX)[i] += w1 * dW;
            ((double*)&this->gradientY)[i] += w2 * dW;
        }
    }
}

float2 Cell::getCenter()
{
	return this->center;
}

float2 Cell::getNode(int i)
{
    return *this->nodes[i];
}

PrimitiveVariable Cell::getPrim()
{
    PrimitiveVariable prim;
    prim.rho = this->cons[0];
    prim.U   = this->cons[1] / this->cons[0];
    prim.V   = this->cons[2] / this->cons[0];

    if( this->interfaceType == incompressible )
        prim.L = 3.0 / 2.0;
    else
        // eq. in GKS Book page 79 at the bottom
        prim.L = (this->fluidParam.K + 2.0)*this->cons[0]
                      / ( 4.0 * ( this->cons[3] - 0.5*(this->cons[1]* this->cons[1] + this->cons[2] * this->cons[2])/this->cons[0] ) );

    return prim;
}

ConservedVariable Cell::getCons()
{
    ConservedVariable tmp;
    tmp.rho  = this->cons[0];
    tmp.rhoU = this->cons[1];
    tmp.rhoV = this->cons[2];
    tmp.rhoE = this->cons[3];
    return tmp;
}

ConservedVariable Cell::getLocalResidual()
{
    return this->residual;
}

ConservedVariable Cell::getUpdate()
{
    return this->updateVal;
}

ConservedVariable Cell::getGradientX()
{
    return this->gradientX;
}

ConservedVariable Cell::getGradientY()
{
    return this->gradientY;
}

Cell * Cell::getNeighborCell(int i)
{
    return this->InterfaceList[i]->getNeigborCell(this);
}

Cell * Cell::getOpposingCell(Interface * askingInterface)
{
    int j = 0;
    for ( int i = 0; i < 4; i++ )
        if ( askingInterface == this->InterfaceList[i] )
            j = i;

    return this->getNeighborCell( (j+2) % 4 );
}

Cell * Cell::findNeighborInDomain()
{
    Cell* neighborCell = NULL;
    // find Interface to domain
    for (int i = 0; i < 4; i++)
    {
        if ( this->InterfaceList[i] != NULL )
            neighborCell = this->InterfaceList[i]->getNeigborCell(this);
    }

    int i = 0;
    return neighborCell;
}

float2 Cell::getConnectivity(int i)
{
    if( this->InterfaceList[i] == NULL) return float2(0.0, 0.0);
    float2 vector;
    vector.x = this->InterfaceList[i]->getCenter().x - this->center.x;
    vector.y = this->InterfaceList[i]->getCenter().y - this->center.y;
    return vector;
}

bool Cell::isGhostCell()
{
    return BoundaryContitionPointer != NULL;
}

double Cell::distance(float2 point)
{
    return sqrt( ( this->center.x - point.x )*( this->center.x - point.x )
               + ( this->center.y - point.y )*( this->center.y - point.y ) );
}

string Cell::toString()
{
	ostringstream tmp;
	tmp << "Center = ( ";
	tmp << this->center.x;
	tmp << " , ";
	tmp << this->center.y;
	tmp << " )" << endl;
    tmp << "volume = " << this->volume << endl;
	return tmp.str();
}

string Cell::valuesToString()
{
    PrimitiveVariable prim = this->getPrim();

    ostringstream tmp;
    tmp << this->toString() << "\n";
    tmp << prim.rho << " " << prim.U << " " << prim.V << " " << prim.L << "\n";
    return tmp.str();
}

string Cell::writeNodes()
{
	ostringstream tmp;
    for(int i = 0; i < 4; i++)
    {
	    tmp << this->nodes[i]->x << " ";
	    tmp << this->nodes[i]->y << " 0.0\n";
    }
	return tmp.str();
}
