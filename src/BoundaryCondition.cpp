#include "BoundaryCondition.h"
#include "GKSSolver.h"
#include "mshReader.h"
#include <sstream>

void BoundaryCondition::addCell(idType id)
{
    this->Cell.push_back(id);
    this->Cell.shrink_to_fit();
}

void BoundaryCondition::findNeighborCells(mshReader& reader)
{
    this->NeighborCell.resize( this->Cell.size() );
    for( int cell = 0; cell < this->Cell.size(); ++cell )
    {
        idType face = reader.Cell2Face[ this->Cell[cell] ][0];

        if( reader.Cell2BC[ reader.Face2Cell[ face ][0] ] == -1 )
            this->NeighborCell[cell] = reader.Face2Cell[ face ][0];
        else
            this->NeighborCell[cell] = reader.Face2Cell[ face ][1];
    }
}

void BoundaryCondition::setGhostCells(GKSSolver & solver)
{    
    #pragma omp parallel for
    for( int cell = 0; cell < Cell.size(); ++cell )
    {
        this->setGhostCell(solver, cell);
    }
}

// ============================================================================

bcWall::bcWall(double U, double V) : U(U), V(V) {}

bcIsothermalWall::bcIsothermalWall(double U, double V, double L) : U(U), V(V), L(L) {}

bcPeriodicGhost::bcPeriodicGhost() {}

void bcPeriodicGhost::findNeighborCells(mshReader & reader)
{
    this->NeighborCell.resize( this->Cell.size() );
    for( int cell = 0; cell < this->Cell.size(); ++cell )
    {
        idType face = reader.Cell2Face[ this->Cell[cell] ][0];
        idType periodicFace = reader.findPeriodicInterface( face );

        if( reader.Cell2BC[ reader.Face2Cell[periodicFace][0] ] == -1)
            this->NeighborCell[cell] = reader.Face2Cell[periodicFace][0];
        else
            this->NeighborCell[cell] = reader.Face2Cell[periodicFace][1];
    }
}

bcInflowParabolic::bcInflowParabolic(double U, double V, double L, Vec2 p0, Vec2 p1)
                                    : U(U), V(V), L(L), p0(p0), p1(p1) {}

bcInflowUniform::bcInflowUniform(double U, double V, double L) : U(U), V(V), L(L) {}

bcOutflow::bcOutflow() {}

// ============================================================================

void bcWall::setGhostCell(GKSSolver & solver, idType cell)
{
    PrimitiveVariable primNeighbor = solver.getPrim( NeighborCell[cell] );
    PrimitiveVariable prim;

    prim.rho = primNeighbor.rho;
    prim.U   = 2.0*this->U - primNeighbor.U;
    prim.V   = 2.0*this->V - primNeighbor.V;
    prim.L   = primNeighbor.L;
        
    solver.setData(Cell[cell], prim);
}

void bcIsothermalWall::setGhostCell(GKSSolver & solver, idType cell)
{
    PrimitiveVariable primNeighbor = solver.getPrim( NeighborCell[cell] );
    PrimitiveVariable prim;

    prim.rho = primNeighbor.rho * ( 2.0*this->L - primNeighbor.L ) / primNeighbor.L;
    prim.U   = 2.0*this->U - primNeighbor.U;
    prim.V   = 2.0*this->V - primNeighbor.V;
    prim.L   = 2.0*this->L - primNeighbor.L;
        
    solver.setData(Cell[cell], prim);
}

void bcPeriodicGhost::setGhostCell(GKSSolver & solver, idType cell)
{
    PrimitiveVariable primNeighbor = solver.getPrim( NeighborCell[cell] );
    PrimitiveVariable prim;
    
    prim.rho = primNeighbor.rho;
    prim.U   = primNeighbor.U;
    prim.V   = primNeighbor.V;
    prim.L   = primNeighbor.L;
        
    solver.setData(Cell[cell], prim);
}

void bcInflowParabolic::setGhostCell(GKSSolver & solver, idType cell)
{
    PrimitiveVariable primNeighbor = solver.getPrim( NeighborCell[cell] );
    PrimitiveVariable prim;

    Vec2 bcVector( this->p1.x - p0.x, this->p1.y - p0.y );

    double H = sqrt( bcVector.x * bcVector.x + bcVector.y * bcVector.y );

    Vec2 cellCenter = solver.getCellCenter( this->NeighborCell[cell] );
    
    Vec2 cellVector( cellCenter.x - this->p0.x, cellCenter.y - this->p0.y );

    double z = ( cellVector.x * bcVector.x + cellVector.y * bcVector.y ) / ( H * H );

    double U_Inflow = 4.0 * this->U * z * ( 1.0 - z );
    double V_Inflow = 4.0 * this->V * z * ( 1.0 - z );

    prim.rho = primNeighbor.rho * ( 2.0*this->L - primNeighbor.L ) / primNeighbor.L;
    prim.U   = 2.0*U_Inflow - primNeighbor.U;
    prim.V   = 2.0*V_Inflow - primNeighbor.V;
    prim.L   = 2.0*this->L  - primNeighbor.L;
        
    solver.setData(Cell[cell], prim);
}

void bcInflowUniform::setGhostCell(GKSSolver & solver, idType cell)
{
    PrimitiveVariable primNeighbor = solver.getPrim( NeighborCell[cell] );
    PrimitiveVariable prim;

    prim.rho = primNeighbor.rho * ( 2.0*this->L - primNeighbor.L ) / primNeighbor.L;
    prim.U   = 2.0*this->U - primNeighbor.U;
    prim.V   = 2.0*this->V - primNeighbor.V;
    prim.L   = 2.0*this->L - primNeighbor.L;
        
    solver.setData(Cell[cell], prim);
}

void bcOutflow::setGhostCell(GKSSolver & solver, idType cell)
{
    PrimitiveVariable primNeighbor = solver.getPrim( NeighborCell[cell] );
    PrimitiveVariable prim;

    prim.rho = primNeighbor.rho;
    prim.U   = primNeighbor.U;
    prim.V   = primNeighbor.V;
    prim.L   = primNeighbor.L;
        
    solver.setData(Cell[cell], prim);
}

// ============================================================================

void bcWall::setGradientGhostCells(GKSSolver & solver) {return;}

void bcIsothermalWall::setGradientGhostCells(GKSSolver & solver) {return;}

void bcPeriodicGhost::setGradientGhostCells(GKSSolver & solver)
{
    #pragma omp parallel for
    for( int cell = 0; cell < Cell.size(); ++cell )
    {
        solver.setCellGradientX( Cell[cell], solver.getCellGradientX( this->NeighborCell[cell] ) );
        solver.setCellGradientY( Cell[cell], solver.getCellGradientY( this->NeighborCell[cell] ) );
    }
}

void bcInflowParabolic::setGradientGhostCells(GKSSolver & solver) {return;}

void bcInflowUniform::setGradientGhostCells(GKSSolver & solver) {return;}

void bcOutflow::setGradientGhostCells(GKSSolver & solver) {return;}
