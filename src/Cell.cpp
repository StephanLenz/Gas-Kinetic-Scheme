
#include "Cell.h"
#include "Interface.h"
#include <sstream>
#include <iostream>
#include <algorithm>    // min()

using namespace std;

Cell::Cell()
{
	this->InterfaceList = new Interface*[4];
    memset(InterfaceList, NULL, 4 * sizeof(Interface*));
}

Cell::Cell(InterfaceType interfaceType, double centerX, double centerY, double dx, double dy, BoundaryCondition* BC, FluidParameter fluidParam)
{
	this->InterfaceList = new Interface*[4];
    memset(InterfaceList, NULL, 4 * sizeof(Interface*));

	this->centerX = centerX;
	this->centerY = centerY;

	this->dx = dx;
	this->dy = dy;

    this->BoundaryContitionPointer = BC;

    this->fluidParam = fluidParam;

    this->interfaceType = interfaceType;
}


Cell::~Cell()
{
	delete [] InterfaceList;
}

void Cell::update(double dt)
{
    double cons_old[4];
    cons_old[0] = this->cons[0];
    cons_old[1] = this->cons[1];
    cons_old[2] = this->cons[2];
    cons_old[3] = this->cons[3];

    // negative interfaces = in flux
    // positive interfaces = out flux
    this->cons[0] += ( this->InterfaceList[0]->getTimeIntegratedFlux().rho
                     + this->InterfaceList[1]->getTimeIntegratedFlux().rho
                     - this->InterfaceList[2]->getTimeIntegratedFlux().rho
                     - this->InterfaceList[3]->getTimeIntegratedFlux().rho
                     ) / (this->dx*this->dy);

    this->cons[1] += ( this->InterfaceList[0]->getTimeIntegratedFlux().rhoU
                     + this->InterfaceList[1]->getTimeIntegratedFlux().rhoU
                     - this->InterfaceList[2]->getTimeIntegratedFlux().rhoU
                     - this->InterfaceList[3]->getTimeIntegratedFlux().rhoU
                     ) / (this->dx*this->dy);

    this->cons[2] += ( this->InterfaceList[0]->getTimeIntegratedFlux().rhoV
                     + this->InterfaceList[1]->getTimeIntegratedFlux().rhoV
                     - this->InterfaceList[2]->getTimeIntegratedFlux().rhoV
                     - this->InterfaceList[3]->getTimeIntegratedFlux().rhoV
                     ) / (this->dx*this->dy);

    this->cons[3] += ( this->InterfaceList[0]->getTimeIntegratedFlux().rhoE
                     + this->InterfaceList[1]->getTimeIntegratedFlux().rhoE
                     - this->InterfaceList[2]->getTimeIntegratedFlux().rhoE
                     - this->InterfaceList[3]->getTimeIntegratedFlux().rhoE
                     ) / (this->dx*this->dy);

    //// ---------------------- FIX ---------------------------------------------
    //if( fabs(this->cons[2]) < 1.0e-12  ) this->cons[2] = 0.0;
    //// ---------------------- FIX ---------------------------------------------

    // Apply Forcing
    this->cons[1] += dt * cons_old[0] * this->fluidParam.Force.x;
    this->cons[2] += dt * cons_old[0] * this->fluidParam.Force.y;

    // compute primary Variables
    this->computePrim();

    // Set the residual for zero values to zero
    //this->residual.rho  = (fabs(cons_old[0]) > 1.0e-12) ? fabs(this->cons[0] - cons_old[0]) / fabs(cons_old[0]) : 0.0;
    //this->residual.rhoU = (fabs(cons_old[1]) > 1.0e-12) ? fabs(this->cons[1] - cons_old[1]) / fabs(cons_old[1]) : 0.0;
    //this->residual.rhoV = (fabs(cons_old[2]) > 1.0e-12) ? fabs(this->cons[2] - cons_old[2]) / fabs(cons_old[2]) : 0.0;
    //this->residual.rhoE = (fabs(cons_old[3]) > 1.0e-12) ? fabs(this->cons[3] - cons_old[3]) / fabs(cons_old[3]) : 0.0;

    this->residual.rho  = fabs(this->cons[0] - cons_old[0]);
    this->residual.rhoU = fabs(this->cons[1] - cons_old[1]);
    this->residual.rhoV = fabs(this->cons[2] - cons_old[2]);
    this->residual.rhoE = fabs(this->cons[3] - cons_old[3]);

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
            double localValue = 4.0 * value * ( this->centerY - this->centerY*this->centerY );
            this->prim[i] = 2.0*localValue - neighborCell->prim[i];
        }
    }
    
    this->computeCons();
}

void Cell::addInterface(Interface* newInterface, int direction)
{
	this->InterfaceList[direction] = newInterface;
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
    //                  CFL-Condition as in Weidong Li's Code
    // ========================================================================
    // The speed of sound is fixed as 1/sqrt(3) here
    // The use of velocity squares is unclear
    // The inclusion of the viscosity is not really complete (compare to code below)

    //double velocitySquare = this->getPrim().U*this->getPrim().U
    //                      + this->getPrim().V*this->getPrim().V;
    //double localTimestep =  min(dx, dy) 
    //                     / ( velocitySquare
    //                       + 1.0/sqrt(3.0) 
    //                       + 2.0*this->fluidParam.nu/min(dx, dy) );


    // ========================================================================

    //double localTimestep = min(dx, dy) / ( max( fabs(this->getPrim().U), fabs(this->getPrim().V) ) 
    //                                     + sqrt( 5.0/(3.0*2.0*this->getPrim().L) )                  // c_s = sqrt( kappa RT ) = sqrt( 5/3 * 1/2lambda )
    //                                     + 2.0*this->fluidParam.nu/min(dx, dy)
    //                                     );

    // ========================================================================
    //                  CFL-Condition as in Guo, Liu et al. (2008)
    // ========================================================================
    // The formular for the speed of sound in Guo, Liu et al (2008) differs
    // from the one at Wikipedia sqrt(RT) != sqrt(kappa RT)

    //double U_max = max( fabs(this->getPrim().U), fabs(this->getPrim().V) );
    double U_max = sqrt(this->getPrim().U*this->getPrim().U + this->getPrim().V*this->getPrim().V);
    double c_s   = sqrt( 1.0 / ( 2.0*this->getPrim().L ) );                      // c_s = sqrt(RT) = c_s = sqrt(1/2lambda)
    double Re    = U_max * min(dx, dy) / this->fluidParam.nu;

    //double localTimestep = min(dx, dy) / ( ( U_max + c_s ) * ( 1.0 + 2.0 / Re ) );
    double localTimestep = min(dx, dy) / ( ( U_max + c_s + 2.0*this->fluidParam.nu / min(dx, dy) ) );

    // ========================================================================
    
    return localTimestep;
}

float2 Cell::getCenter()
{
	float2 center;
	center.x = this->centerX;
	center.y = this->centerY;
	return center;
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

float2 Cell::getDx()
{
    float2 tmp;
    tmp.x = this->dx;
    tmp.y = this->dy;
    return tmp;
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
        if ((this->InterfaceList[i] != NULL) && !this->InterfaceList[i]->getNeigborCell(this)->isGhostCell())
            neighborCell = this->InterfaceList[i]->getNeigborCell(this);
    }
    return neighborCell;
}

bool Cell::isGhostCell()
{
    return !(BoundaryContitionPointer== NULL);
}

string Cell::toString()
{
	ostringstream tmp;
	tmp << "Center = ( ";
	tmp << this->centerX;
	tmp << " , ";
	tmp << this->centerY;
	tmp << " )";
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
	tmp << this->centerX - 0.5*dx << " ";
	tmp << this->centerY - 0.5*dy << " 0.0\n";

	tmp << this->centerX + 0.5*dx << " ";
	tmp << this->centerY - 0.5*dy << " 0.0\n";

	tmp << this->centerX + 0.5*dx << " ";
	tmp << this->centerY + 0.5*dy << " 0.0\n";

	tmp << this->centerX - 0.5*dx << " ";
	tmp << this->centerY + 0.5*dy << " 0.0\n";
	return tmp.str();
}
