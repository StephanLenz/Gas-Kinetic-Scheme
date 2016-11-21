#include "BoundaryCondition.h"
#include "GKSSolver.h"
#include <sstream>

BoundaryCondition::BoundaryCondition()
{
}

BoundaryCondition::BoundaryCondition( BoundaryConditionType type,
                                      double rho, double U, double V, double L)
{
    this->type = type;
    this->value.rho = rho;
    this->value.U = U;
    this->value.V = V;
    this->value.L = L;
}

BoundaryCondition::~BoundaryCondition()
{
}

BoundaryConditionType BoundaryCondition::getType()
{
    return this->type;
}

PrimitiveVariable BoundaryCondition::getValue()
{
    return value;
}

void BoundaryCondition::addCell(idType id)
{
    this->Cell.push_back(id);
    this->Cell.shrink_to_fit();
}

void BoundaryCondition::addNeighborCell(idType id)
{
    this->NeighborCell.push_back(id);
    this->NeighborCell.shrink_to_fit();
}

void BoundaryCondition::setGhostCells(GKSSolver & solver)
{
    #pragma omp parallel for
    for( int i = 0; i < Cell.size(); ++i )
    {

        PrimitiveVariable primNeighbor = solver.getPrim( NeighborCell[i] );
        PrimitiveVariable prim;

        switch ( this->type )
        {
            // ====================================================================
            case wall:
            {
                prim.rho = primNeighbor.rho;
                prim.U   = 2.0*this->value.U - primNeighbor.U;
                prim.V   = 2.0*this->value.V - primNeighbor.V;
                prim.L   = primNeighbor.L;

                break;
            }
            // ====================================================================
            case isothermalWall:
            {
                prim.rho = primNeighbor.rho * ( 2.0*this->value.L - primNeighbor.L ) / primNeighbor.L;
                prim.U   = 2.0*this->value.U - primNeighbor.U;
                prim.V   = 2.0*this->value.V - primNeighbor.V;
                prim.L   = 2.0*this->value.L - primNeighbor.L;
                
                break;
            }
            // ====================================================================
            case inlet:
            {
                prim.rho = primNeighbor.rho;
                prim.U   = 2.0*this->value.U - primNeighbor.U;
                prim.V   = 2.0*this->value.V - primNeighbor.V;
                prim.L   = primNeighbor.L;

                break;
            }
            // ====================================================================
            case outlet:
            {
                prim.rho = primNeighbor.rho;
                prim.U   = primNeighbor.U;
                prim.V   = primNeighbor.V;
                prim.L   = primNeighbor.L;

                break;
            }
            // ====================================================================
        }
        
        solver.setData(Cell[i], prim);
    }
}

std::string BoundaryCondition::toString()
{
    std::stringstream tmp;
    if      ( this->type == wall )
        tmp << "Wall, U = " << this->value.U << ", V = " << this->value.V;
    else if ( this->type == isothermalWall )
        tmp << "Isothermal Wall, U = " << this->value.U << ", V = " << this->value.V << ", lambda = " << this->value.L;
    else if ( this->type == periodic )
        tmp << "Periodic Boundary";
    return tmp.str();
}
