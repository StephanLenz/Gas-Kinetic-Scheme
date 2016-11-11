// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      Cell.h
//
//      Function:
//          Storage of cell data
//          Implementation of cell related computations
//              Update of Conserved Variables
//              Forcing
//              Timestep computations
//
// ============================================================================

#include "Cell.h"
#include "Interface.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>    // min()

using namespace std;

unsigned long int Cell::numberOfCells = 1;

Cell::Cell()
{
    InterfaceList[0] = NULL;
    InterfaceList[1] = NULL;
    InterfaceList[2] = NULL;
    InterfaceList[3] = NULL;
}

// ============================================================================
//      This constructor initializes the Cell.
//      It also computes:
//          Centroid
//          Volume
//
//      Parameters:
//          interfaceType:  compressible / incompressible
//          nodes:          pointer to an array with pointers to the four nodes
//                          of the cell order counter clockwise
//          BC:             Boundary Condition Pointer
//          fluidParam:     FluidParameter Object
// ============================================================================
Cell::Cell(InterfaceType interfaceType, Node** nodes, BoundaryCondition* BC, FluidParameter fluidParam)
    : BoundaryContitionPointer(BC)
{
    // ========================================================================
    //                  Copy attributes
    // ========================================================================
    InterfaceList[0] = NULL;
    InterfaceList[1] = NULL;
    InterfaceList[2] = NULL;
    InterfaceList[3] = NULL;

    this->nInterfaces = 0;

    for(int i = 0; i < 4; i++)
        this->nodes[i] = nodes[i];

    this->fluidParam = fluidParam;

    this->interfaceType = interfaceType;
    // ========================================================================
    //                  Copy attributes
    // ========================================================================
    this->ID = Cell::numberOfCells++;

    // ========================================================================
    
    // ========================================================================
    //                  Compute Centroid and Volume of the Quad
    // ========================================================================
    Node triCenter[2];
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
    // ========================================================================

    //// ========================================================================
    ////                  Compute Volume of the Quad (old)
    //// ========================================================================
    //this->volume = 0.5 * fabs( this->nodes[0]->x * ( this->nodes[1]->y - this->nodes[3]->y ) 
    //                         + this->nodes[1]->x * ( this->nodes[3]->y - this->nodes[0]->y ) 
    //                         + this->nodes[3]->x * ( this->nodes[0]->y - this->nodes[1]->y ) )
    //             + 0.5 * fabs( this->nodes[2]->x * ( this->nodes[3]->y - this->nodes[1]->y ) 
    //                         + this->nodes[3]->x * ( this->nodes[1]->y - this->nodes[2]->y ) 
    //                         + this->nodes[1]->x * ( this->nodes[2]->y - this->nodes[3]->y ) );
    //// ========================================================================
    //
    //// ========================================================================
    ////                  Compute geometric Center of the Quad (Old)
    //// ========================================================================
    //this->center.x = 0.25 * (this->nodes[0]->x + this->nodes[1]->x + this->nodes[2]->x + this->nodes[3]->x);
    //this->center.y = 0.25 * (this->nodes[0]->y + this->nodes[1]->y + this->nodes[2]->y + this->nodes[3]->y);
    //// ========================================================================
    
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
//      This method returns the local timestep computed by an CFL like condition
// ============================================================================
double Cell::getLocalTimestep()
{
    PrimitiveVariable prim = this->getPrim();

    double U_max = sqrt(prim.U*prim.U + prim.V*prim.V);
    double c_s   = sqrt( 1.0 / ( 2.0*prim.L ) );           // c_s = sqrt(RT) = c_s = sqrt(1/2lambda)

    double localTimestep = minDx / ( U_max + c_s + 2.0*this->fluidParam.nu / minDx );
    
    return localTimestep;
}

// ============================================================================
//      This method applies the forcing to the cell and therefore increases
//      the momentum.
//
//      Parameters:
//          dt:   global time step
// ============================================================================
void Cell::applyForcing(double dt)
{
    Node Force;
    
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

    int breakpoint = 0;
    // ========================================================================

    // ========================================================================
    //                  compute new total Energy with increased momentum
    // ========================================================================
    this->cons[3] = prim.rho * (this->fluidParam.K + 2.0) / (4.0*prim.L)
                  + 0.5 * ( this->cons[1] * this->cons[1] 
                          + this->cons[2] * this->cons[2]
                          ) / prim.rho;
    // ========================================================================
}

// ============================================================================
//      This method computes the primitive varables in the ghost cells 
//      according to the boundary conditions.
//      Afterwards the primitive variables are transformed to conserved
//      variables.
// ============================================================================
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
        // ====================================================================
        case periodicGhost:
        {
            prim = this->InterfaceList[1]->getNeigborCell(this)->getPrim();
            this->gradientX = this->InterfaceList[1]->getNeigborCell(this)->getGradientX();
            this->gradientY = this->InterfaceList[1]->getNeigborCell(this)->getGradientY();
        }
        // ====================================================================
    }
    
    this->computeCons(prim);
}

// ============================================================================
//      This method computes the conserved variables in this cell from the 
//      struct of primitive variables given as a parameter.
//
//      Parameters:
//          prim:   struct of primitive Variables
// ============================================================================
void Cell::computeCons(PrimitiveVariable prim)
{
    this->cons[0] = prim.rho;
    this->cons[1] = prim.rho * prim.U; 
    this->cons[2] = prim.rho * prim.V;
    this->cons[3] = prim.rho * (this->fluidParam.K + 2.0) / (4.0*prim.L)
                  + 0.5 * ( this->cons[1] * this->cons[1] 
                          + this->cons[2] * this->cons[2]
                          )/prim.rho;
}

// ============================================================================
//      This method computes the gradients of the conserved variables per
//      Cell by the Least square Methods as described in Blazek (Comp. Fluid
//      Dynamics: Principles and Applications)
//      
//      Since the Ghost cells have only one neighbor, the Gradient is computed 
//      by finite difference as a direction derivative
// ============================================================================
void Cell::computeLeastSquareGradients()
{
    if( this->isGhostCell() )
    {
        // loop over conserved variables
        for(int i = 0; i < 4; i++)
        {
            double dx = this->InterfaceList[0]->getNeigborCell(this)->getCenter().x - this->center.x;
            double dy = this->InterfaceList[0]->getNeigborCell(this)->getCenter().y - this->center.y;

            double dW = ((double*)&this->InterfaceList[0]->getNeigborCell(this)->getCons())[i] - this->cons[i];

            ((double*)&this->gradientX)[i] = dW / (dx*dx + dy*dy) * dx;
            ((double*)&this->gradientY)[i] = dW / (dx*dx + dy*dy) * dy;
        }
        return;
    }

    // loop over conserved variables
    for(int i = 0; i < 4; i++)
    {
        ((double*)&this->gradientX)[i] = 0.0;
        ((double*)&this->gradientY)[i] = 0.0;

        // loop over 4 adjacent cells
        for(int j = 0; j < 4; j++)
        {
            double dx = this->InterfaceList[j]->getNeigborCell(this)->getCenter().x - this->center.x;
            double dy = this->InterfaceList[j]->getNeigborCell(this)->getCenter().y - this->center.y;

            // distance weighting
            double weight = 1.0 / sqrt( dx*dx + dy*dy );
            dx *= weight;
            dy *= weight;

            double dW = weight * ( ((double*)&this->InterfaceList[j]->getNeigborCell(this)->getCons())[i] - this->cons[i] );

            double w1 = dx / (r11*r11) - r12/r11 * ( dy - r12/r11 * dx ) / (r22*r22);
            double w2 =                            ( dy - r12/r11 * dx ) / (r22*r22);

            ((double*)&this->gradientX)[i] += w1 * dW;
            ((double*)&this->gradientY)[i] += w2 * dW;
        }
    }
}

// ============================================================================
//      This method updates the conserved variables in the cell.
//      It also computes the absolute residual change in the cell.
//
//      Parameters:
//          dt:   global time step
// ============================================================================
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
    //                        Compute Residual and store new values
    // ========================================================================
    this->residual.rho  = fabs(this->cons[0] - this->cons_old[0]);
    this->residual.rhoU = fabs(this->cons[1] - this->cons_old[1]);
    this->residual.rhoV = fabs(this->cons[2] - this->cons_old[2]);
    this->residual.rhoE = fabs(this->cons[3] - this->cons_old[3]);

    if( this->ID == 1 )
        int breakPoint = 0;

    // store values of this time step for residual computation in the next one
    this->cons_old[0] = this->cons[0];
    this->cons_old[1] = this->cons[1];
    this->cons_old[2] = this->cons[2];
    this->cons_old[3] = this->cons[3];
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

    // ========================================================================

    ConservedVariable update_1;

    update_1.rho     = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_1().rho
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_1().rho
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_1().rho
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_1().rho
                     ) / this->volume;

    update_1.rhoU    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_1().rhoU
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_1().rhoU
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_1().rhoU
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_1().rhoU
                     ) / this->volume;

    update_1.rhoV    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_1().rhoV
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_1().rhoV
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_1().rhoV
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_1().rhoV
                     ) / this->volume;

    update_1.rhoE    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_1().rhoE
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_1().rhoE
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_1().rhoE
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_1().rhoE
                     ) / this->volume;

    ConservedVariable update_2;

    update_2.rho     = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_2().rho
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_2().rho
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_2().rho
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_2().rho
                     ) / this->volume;

    update_2.rhoU    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_2().rhoU
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_2().rhoU
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_2().rhoU
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_2().rhoU
                     ) / this->volume;

    update_2.rhoV    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_2().rhoV
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_2().rhoV
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_2().rhoV
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_2().rhoV
                     ) / this->volume;

    update_2.rhoE    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_2().rhoE
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_2().rhoE
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_2().rhoE
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_2().rhoE
                     ) / this->volume;

    ConservedVariable update_3;

    update_3.rho     = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_3().rho
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_3().rho
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_3().rho
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_3().rho
                     ) / this->volume;

    update_3.rhoU    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_3().rhoU
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_3().rhoU
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_3().rhoU
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_3().rhoU
                     ) / this->volume;

    update_3.rhoV    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_3().rhoV
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_3().rhoV
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_3().rhoV
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_3().rhoV
                     ) / this->volume;

    update_3.rhoE    = ( Sign0 * this->InterfaceList[0]->getTimeIntegratedFlux_3().rhoE
                       + Sign1 * this->InterfaceList[1]->getTimeIntegratedFlux_3().rhoE
                       + Sign2 * this->InterfaceList[2]->getTimeIntegratedFlux_3().rhoE
                       + Sign3 * this->InterfaceList[3]->getTimeIntegratedFlux_3().rhoE
                     ) / this->volume;

    // ========================================================================

    this->updateVal.rho  = update.rho ;
    this->updateVal.rhoU = update.rhoU;
    this->updateVal.rhoV = update.rhoV;
    this->updateVal.rhoE = update.rhoE;

    // ========================================================================

    // ========================================================================
    //                        Output Cell updates for smaller cell
    // ========================================================================
    //ofstream file;
    //file.precision(16);
    //file.open("out/update_rhoV.dat", fstream::app );
    //file << scientific << setw(25) << setfill(' ') << update_1.rhoV;
    //file << scientific << setw(25) << setfill(' ') << update_2.rhoV;
    //file << scientific << setw(25) << setfill(' ') << update_3.rhoV << endl;
    // ========================================================================
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          initialization Methods
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method add one interface to the cell and increases the interface
//      counter for this cell.
//
//      Parameters:
//          newInterface:   Pointer to the Interface to be added
// ============================================================================
void Cell::addInterface(Interface* newInterface)
{
	this->InterfaceList[this->nInterfaces++] = newInterface;
}

// ============================================================================
//      This method computes the minimal distance in the cell. this value is
//      needed for the timestep computation.
// ============================================================================
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
                    Node node( this->nodes[j]->x, this->nodes[j]->y );
                    this->minDx = min( this->minDx, this->InterfaceList[i]->distance(node) );
                }
            }
        }
    }

    int breakPoint = 0;
}

// ============================================================================
//      This method computes the least square coefficients for this cell
//      according to Blazek (Comp. Fluid Dynamics: Principles and Applications)
// ============================================================================
void Cell::computeLeastSquareCoefficients()
{
    this->r11 = 0.0;
    this->r12 = 0.0;
    this->r22 = 0.0;

    // loop over all adjacent Cells
    for(int i = 0; i < 4; i++)
    {
        double dx = this->InterfaceList[i]->getNeigborCell(this)->getCenter().x - this->center.x;
        double dy = this->InterfaceList[i]->getNeigborCell(this)->getCenter().y - this->center.y;

        // distance weighting
        double weight = 1.0 / sqrt( dx*dx + dy*dy );
        dx *= weight;
        dy *= weight;

        this->r11 += dx*dx;
        this->r12 += dx*dy;
        this->r22 += dy*dy;
    }

    this->r11 = sqrt(this->r11);
    this->r12 = this->r12 / r11;
    this->r22 = sqrt(this->r22 - r12*r12);
}

// ============================================================================
//      This method computes the conserved variables according to the
//      primitive Variables passed as parameters.
//      It also assignes the old conserved variables
//
//      Parameters:
//          rho:        Density
//          U:          x-velocity
//          V:          y-velocity
//          L:          Lambda = 1/2RT
// ============================================================================
void Cell::setValues(double rho, double U, double V, double L)
{
    PrimitiveVariable prim;
	prim.rho = rho;
	prim.U   = U;
	prim.V   = V;
	prim.L   = L;

    this->computeCons(prim);

    this->cons_old[0] = this->cons[0];
    this->cons_old[1] = this->cons[1];
    this->cons_old[2] = this->cons[2];
    this->cons_old[3] = this->cons[3];
}

void Cell::setCons(ConservedVariable cons)
{
    this->cons[0] = cons.rho;
    this->cons[1] = cons.rhoU;
    this->cons[2] = cons.rhoV;
    this->cons[3] = cons.rhoE;
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

// ============================================================================
//      This method returns the Cell ID.
// ============================================================================
unsigned long int Cell::getID()
{
    return this->ID;
}

// ============================================================================
//      This method returns the location of the i-th node of the cell.
//
//      Parameters:
//          i:   the local id of the node to be returned
// ============================================================================
Node Cell::getNode(int i)
{
    return *this->nodes[i];
}

Interface * Cell::getInterface(int i)
{
    return this->InterfaceList[i];
}

// ============================================================================
//      This method returns a pointer to the adjacent cell connected by the 
//      i-th interface
//
//      Parameters:
//          i:   The local id of the interface that connects the returned cell.
// ============================================================================
Cell * Cell::getNeighborCell(int i)
{
    return this->InterfaceList[i]->getNeigborCell(this);
}

// ============================================================================
//      This method is designed for ghost cells. In case of a ghost cell it
//      it returns the corresponding cell in the domain.
// ============================================================================
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

// ============================================================================
//      This method returns a vector from the cell center to the center of the
//      i-th interface.
// ============================================================================
Node Cell::getConnectivity(int i)
{
    if( this->InterfaceList[i] == NULL) return Node(0.0, 0.0);
    Node vector;
    vector.x = this->InterfaceList[i]->getCenter().x - this->center.x;
    vector.y = this->InterfaceList[i]->getCenter().y - this->center.y;
    return vector;
}

// ============================================================================
//      This method returns wether the cell is a ghost cells or not.
//      This decision is done on basis ofthe boundary condition pointer.
//      Only shost cells have a pointer to a boundary condition.
// ============================================================================
bool Cell::isGhostCell()
{
    return BoundaryContitionPointer != NULL;
}

BoundaryCondition * Cell::getBoundaryConditionPointer()
{
    return this->BoundaryContitionPointer;
}

// ====================================================================================================================
//
//                          get Geometry
//
// ====================================================================================================================

// ============================================================================
//      This method returns the cell volume.
// ============================================================================
double Cell::getVolume()
{
    return this->volume;
}

// ============================================================================
//      This method returns the cell center.
// ============================================================================
Node Cell::getCenter()
{
	return this->center;
}

double Cell::getMinDx()
{
    return this->minDx;
}

// ============================================================================
//      This method returns the distance from the cell center of this cell
//      to the point passed as parameter.
//
//      Parameters:
//          point:   Point to which the distance is computed
// ============================================================================
double Cell::distance(Node point)
{
    return sqrt( ( this->center.x - point.x )*( this->center.x - point.x )
               + ( this->center.y - point.y )*( this->center.y - point.y ) );
}

// ====================================================================================================================
//
//                          get Data
//
// ====================================================================================================================

// ============================================================================
//      This method returns the primitve variabeles in this cell.
//      In case of incompressible Interfaces lambda = 1.5 is returned.
// ============================================================================
PrimitiveVariable Cell::getPrim()
{
    PrimitiveVariable prim;
    prim.rho = this->cons[0];
    prim.U   = this->cons[1] / this->cons[0];
    prim.V   = this->cons[2] / this->cons[0];

    if( this->interfaceType == incompressible )
        prim.L = 3.0 / 2.0;
    else
        prim.L = (this->fluidParam.K + 2.0)*this->cons[0]
                      / ( 4.0 * ( this->cons[3] - 0.5*(this->cons[1] * this->cons[1] 
                                                     + this->cons[2] * this->cons[2])/this->cons[0] ) );

    return prim;
}

// ============================================================================
//      This method returns the conserved variables in this cell
// ============================================================================
ConservedVariable Cell::getCons()
{
    ConservedVariable tmp;
    tmp.rho  = this->cons[0];
    tmp.rhoU = this->cons[1];
    tmp.rhoV = this->cons[2];
    tmp.rhoE = this->cons[3];
    return tmp;
}

ConservedVariable Cell::getConsOld()
{
    ConservedVariable tmp;
    tmp.rho  = this->cons_old[0];
    tmp.rhoU = this->cons_old[1];
    tmp.rhoV = this->cons_old[2];
    tmp.rhoE = this->cons_old[3];
    return tmp;
}

// ============================================================================
//      This method returns the local absolute residual change.
// ============================================================================
ConservedVariable Cell::getLocalResidual()
{
    return this->residual;
}

// ============================================================================
//      This method returns the update of conserved variables in this time step.
// ============================================================================
ConservedVariable Cell::getUpdate()
{
    return this->updateVal;
}

// ============================================================================
//      This method returns the X-Component of the Gradients of the conserved
//      Variables.
// ============================================================================
ConservedVariable Cell::getGradientX()
{
    return this->gradientX;
}

// ============================================================================
//      This method returns the Y-Component of the Gradients of the conserved
//      Variables.
// ============================================================================
ConservedVariable Cell::getGradientY()
{
    return this->gradientY;
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
//      This method returns a short string with information on the Cell.
// ============================================================================
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

// ============================================================================
//      This method returns the a string with short information on the cell
//      and includes the primitive variables in this cell.
// ============================================================================
string Cell::valuesToString()
{
    PrimitiveVariable prim = this->getPrim();

    ostringstream tmp;
    tmp << this->toString() << "\n";
    tmp << prim.rho << " " << prim.U << " " << prim.V << " " << prim.L << "\n";
    return tmp.str();
}

// ============================================================================
//      This method returns a string with the coordinates of the nodes 
//      defining the cell.
// ============================================================================
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

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================