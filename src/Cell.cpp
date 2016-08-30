
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
    //                  Compute Volume of the Quad
    // ========================================================================
    volume = 0.5 * fabs( this->nodes[0]->x * ( this->nodes[1]->y - this->nodes[3]->y ) 
                      + this->nodes[1]->x * ( this->nodes[3]->y - this->nodes[0]->y ) 
                      + this->nodes[3]->x * ( this->nodes[0]->y - this->nodes[1]->y ) )
           + 0.5 * fabs( this->nodes[2]->x * ( this->nodes[3]->y - this->nodes[1]->y ) 
                      + this->nodes[3]->x * ( this->nodes[1]->y - this->nodes[2]->y ) 
                      + this->nodes[1]->x * ( this->nodes[2]->y - this->nodes[3]->y ) );
    // ========================================================================
    
    // ========================================================================
    //                  Compute geometric Center of the Quad
    // ========================================================================
    this->center.x = 0.25 * (this->nodes[0]->x + this->nodes[1]->x + this->nodes[2]->x + this->nodes[3]->x);
    this->center.y = 0.25 * (this->nodes[0]->y + this->nodes[1]->y + this->nodes[2]->y + this->nodes[3]->y);
    // ========================================================================


    //// ========================================================================
    ////                  Compute maximal length of the Quad
    //// ========================================================================
    //double xMax = -1.0e99;
    //double xMin =  1.0e99;
    //double yMax = -1.0e99;
    //double yMin =  1.0e99;
    //for (int i = 0; i < 4; i++)
    //{
    //    if(this->nodes[i]->x > xMax)    xMax = this->nodes[i]->x;
    //    if(this->nodes[i]->x < xMin)    xMin = this->nodes[i]->x;
    //    if(this->nodes[i]->y > yMax)    yMax = this->nodes[i]->y;
    //    if(this->nodes[i]->y < yMin)    yMin = this->nodes[i]->y;
    //}
    //this->dx.x = xMax - xMin;
    //this->dx.y = yMax - yMin;
    //// ========================================================================
    
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
        for(int i = 0; i < 4; i++)
        {
            for( int j = 0; j < 4; j++)
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
    PrimitiveVariable prim = getPrim();

    double U_max = sqrt(prim.U*prim.U + prim.V*prim.V);
    double c_s   = sqrt( 1.0 / ( 2.0*prim.L ) );                      // c_s = sqrt(RT) = c_s = sqrt(1/2lambda)

    // TODO: The minimal distance in a cell is as the minimum of the maximal length is the coordinate directions
    double localTimestep = minDx / ( U_max + c_s + 2.0*this->fluidParam.nu / minDx );

    // ========================================================================
    
    return localTimestep;
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
