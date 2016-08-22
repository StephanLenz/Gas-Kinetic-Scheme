
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
    memset(InterfaceList, NULL, 4 * sizeof(Interface*));

    this->nInterfaces = 0;

    for(int i = 0; i < 4; i++)
        this->nodes[i] = nodes[i];

    volume = 0.5 * abs( this->nodes[0]->x * ( this->nodes[1]->y - this->nodes[3]->y ) 
                      + this->nodes[1]->x * ( this->nodes[3]->y - this->nodes[0]->y ) 
                      + this->nodes[3]->x * ( this->nodes[0]->y - this->nodes[1]->y ) )
           + 0.5 * abs( this->nodes[2]->x * ( this->nodes[3]->y - this->nodes[1]->y ) 
                      + this->nodes[3]->x * ( this->nodes[1]->y - this->nodes[2]->y ) 
                      + this->nodes[1]->x * ( this->nodes[2]->y - this->nodes[3]->y ) );

    this->center.x = 0.25 * (this->nodes[0]->x + this->nodes[1]->x + this->nodes[2]->x + this->nodes[3]->x);
    this->center.y = 0.25 * (this->nodes[0]->y + this->nodes[1]->y + this->nodes[2]->y + this->nodes[3]->y);


    // compute maximal legth of the cell
    double xMax = -1.0e99;
    double xMin =  1.0e99;
    double yMax = -1.0e99;
    double yMin =  1.0e99;
    for (int i = 0; i < 4; i++)
    {
        if(this->nodes[i]->x > xMax)    xMax = this->nodes[i]->x;
        if(this->nodes[i]->x < xMin)    xMin = this->nodes[i]->x;
        if(this->nodes[i]->y > yMax)    yMax = this->nodes[i]->y;
        if(this->nodes[i]->y < yMin)    yMin = this->nodes[i]->y;
    }
    this->dx.x = xMax - xMin;
    this->dx.y = yMax - yMin;


    this->BoundaryContitionPointer = BC;

    this->fluidParam = fluidParam;

    this->interfaceType = interfaceType;
}


Cell::~Cell()
{
}

void Cell::update(double dt)
{

    // ========================================================================
    // negative interfaces = in flux
    // positive interfaces = out flux
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

    ConservedVariable lefFlux = this->InterfaceList[0]->getTimeIntegratedFlux();
    ConservedVariable botFlux = this->InterfaceList[1]->getTimeIntegratedFlux();
    ConservedVariable rigFlux = this->InterfaceList[2]->getTimeIntegratedFlux();
    ConservedVariable topFlux = this->InterfaceList[3]->getTimeIntegratedFlux();

    int i = 0;
    // ========================================================================

    // compute primary Variables
    this->computePrim();

    this->residual.rho  = fabs(this->cons[0] - this->cons_old[0]);
    this->residual.rhoU = fabs(this->cons[1] - this->cons_old[1]);
    this->residual.rhoV = fabs(this->cons[2] - this->cons_old[2]);
    this->residual.rhoE = fabs(this->cons[3] - this->cons_old[3]);

    // store values of this time step for residual computation in the next one
    this->cons_old[0] = this->cons[0];
    this->cons_old[1] = this->cons[1];
    this->cons_old[2] = this->cons[2];
    this->cons_old[3] = this->cons[3];

   int j = 0;

}

void Cell::applyBoundaryCondition()
{

    Cell* neighborCell = this->findNeighborInDomain();

    // if no neighbor was found
    if (neighborCell == NULL)
    {
        // search any interface and take the neigbor 
        for (int i = 0; i < 4; i++)
        {
            if (this->InterfaceList[i] != NULL)
            {
                neighborCell = InterfaceList[i]->getNeigborCell(this)->findNeighborInDomain();
                break;
            }
        }
    }

    // loop over all primitive variables
    for (int i = 0; i < 4; i++)
    {
        int type =     this->BoundaryContitionPointer->getType(i);
        double value = this->BoundaryContitionPointer->getValue(i);

        if (type == 0)
        {
            this->prim[i] = 2.0*value - neighborCell->prim[i];
        }
        else if (type == 1)
        {
            this->prim[i] = neighborCell->prim[i];
        }
        else if ( type == 2 )
        {
            double localValue = 4.0 * value * ( this->center.y - this->center.y*this->center.y );
            this->prim[i] = 2.0*localValue - neighborCell->prim[i];
        }
        else if ( type == 3 )
        {
            double lambda_0 = this->BoundaryContitionPointer->getValue(3);
            double lambda_1 = neighborCell->prim[3];
            this->prim[i] = neighborCell->prim[i] * ( 2.0*lambda_0 - lambda_1 ) / lambda_1;
        }
    }
    
    this->computeCons();
}

void Cell::applyForcing(double dt)
{
    float2 Force;

    Force.x = this->fluidParam.Force.x + this->fluidParam.BoussinesqForce.x * (  this->cons[0] - this->fluidParam.rhoReference ) / this->fluidParam.rhoReference;
    Force.y = this->fluidParam.Force.y + this->fluidParam.BoussinesqForce.y * (  this->cons[0] - this->fluidParam.rhoReference ) / this->fluidParam.rhoReference;

    // Apply Forcing to momentum components
    this->cons[1] += dt * this->cons[0] * Force.x;
    this->cons[2] += dt * this->cons[0] * Force.y;

    // ========== some Tests ==================================================
    double xForce = dt * this->cons[0] * Force.x;
    double yForce = dt * this->cons[0] * Force.y;

    int i = 0;
    // ========================================================================

    // compute new Energy with increased momentum
    this->cons[3] = this->prim[0] * (this->fluidParam.K + 2.0) / (4.0*this->prim[3])
                  + 0.5 * ( this->cons[1] * this->cons[1] + this->cons[2] * this->cons[2] ) / this->prim[0];

    // Compute primitive Variables
    this->computePrim();

    int j = 0;
}

void Cell::addInterface(Interface* newInterface)
{
	this->InterfaceList[this->nInterfaces++] = newInterface;
}

void Cell::setValues(double rho, double u, double v, double L)
{
	this->prim[0] = rho;
	this->prim[1] = u;
	this->prim[2] = v;
	this->prim[3] = L;

    this->computeCons();
}

void Cell::computePrim()
{
    this->prim[0] = this->cons[0];
    this->prim[1] = this->cons[1] / this->cons[0];
    this->prim[2] = this->cons[2] / this->cons[0];

    if( this->interfaceType == incompressible )
        this->prim[3] = 3.0 / 2.0;
    else
        // eq. in GKS Book page 79 at the bottom
        this->prim[3] = (this->fluidParam.K + 2.0)*this->cons[0]
                      / ( 4.0 * ( this->cons[3] - 0.5*(this->cons[1]* this->cons[1] + this->cons[2] * this->cons[2])/this->cons[0] ) );
}

void Cell::computeCons()
{
    this->cons[0] = this->prim[0];
    this->cons[1] = this->prim[0] * this->prim[1]; 
    this->cons[2] = this->prim[0] * this->prim[2];
    // inverse of eq. in GKS Book page 79 at the bottom
    this->cons[3] = this->prim[0] * (this->fluidParam.K + 2.0) / (4.0*this->prim[3])
                  + 0.5 * (this->cons[1] * this->cons[1] + this->cons[2] * this->cons[2])/this->prim[0];
}

double Cell::getLocalTimestep()
{
    // ========================================================================
    //                  CFL-Condition as in Guo, Liu et al. (2008)
    // ========================================================================
    // The formular for the speed of sound in Guo, Liu et al (2008) differs
    // from the one at Wikipedia sqrt(RT) != sqrt(kappa RT)

    //double U_max = max( fabs(this->getPrim().U), fabs(this->getPrim().V) );
    double U_max = sqrt(this->getPrim().U*this->getPrim().U + this->getPrim().V*this->getPrim().V);
    double c_s   = sqrt( 1.0 / ( 2.0*this->getPrim().L ) );                      // c_s = sqrt(RT) = c_s = sqrt(1/2lambda)

    // TODO: The minimal distance in a cell is here computed  as if the cell where a square. This is ofcourse wrong!
    double localTimestep = min(this->dx.x, this->dx.y) / ( ( U_max + c_s + 2.0*this->fluidParam.nu / min(this->dx.x, this->dx.y) ) );

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

PrimaryVariable Cell::getPrim()
{
    PrimaryVariable tmp;
    tmp.rho = this->prim[0];
    tmp.U   = this->prim[1];
    tmp.V   = this->prim[2];
    tmp.L   = this->prim[3];
    return tmp;
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
        //if ((this->InterfaceList[i] != NULL) && !this->InterfaceList[i]->getNeigborCell(this)->isGhostCell())
        if ( this->InterfaceList[i] != NULL )
            neighborCell = this->InterfaceList[i]->getNeigborCell(this);
    }

    int i = 0;
    return neighborCell;
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
    ostringstream tmp;
    tmp << this->toString() << "\n";
    tmp << this->prim[0] << " " << this->prim[1] << " " << this->prim[2] << " " << this->prim[3] << "\n";
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
