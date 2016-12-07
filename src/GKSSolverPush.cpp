
#include "GKSSolverPush.h"
#include "BoundaryCondition.h"
#include "Types.h"
#include "mshReader.h"
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

GKSSolverPush::GKSSolverPush() : GKSSolver()
{
}

GKSSolverPush::GKSSolverPush(Parameters param, FluidParameter fluidParam)
         : GKSSolver(param, fluidParam)
{
}

GKSSolverPush::~GKSSolverPush()
{
}

bool GKSSolverPush::readProblem(string filename)
{
    mshReader reader;

    if( ! reader.readProblem(filename) ) return false;
    
    this->numberOfNodes        = reader.Nodes.size();

    this->numberOfCells        = reader.Cell2Node.size();

    this->numberOfInterfaces   = reader.Face2Node.size();

    // ========================================================================
    //              Allocate Memory
    // ========================================================================

    this->NodeCenter.resize( numberOfNodes );

    this->CellData.resize( numberOfCells );
    this->CellDataOld.resize( numberOfCells );

    this->CellUpdate.resize( numberOfCells );

    this->CellGradientX.resize( numberOfCells );
    this->CellGradientY.resize( numberOfCells );

    this->Cell2Node.resize( numberOfCells );
    this->Cell2Interface.resize( numberOfCells );
    this->CellBoundaryCondition.resize( numberOfCells );

    this->CellCenter.resize( numberOfCells );
    this->CellVolume.resize( numberOfCells );
    this->CellMinDx.resize( numberOfCells );
    this->CellLSCoeff.resize( numberOfCells );

    this->InterfaceFlux.resize( numberOfInterfaces );

    this->Interface2Node.resize( numberOfInterfaces );
    this->Interface2Cell.resize( numberOfInterfaces );

    this->InterfaceCenter.resize( numberOfInterfaces );
    this->InterfaceNormal.resize( numberOfInterfaces );
    this->InterfaceArea.resize( numberOfInterfaces );
    this->InterfaceDistance.resize( numberOfInterfaces );
    this->Interface2CellCenterDistance.resize( numberOfInterfaces );

    this->InterfaceAdd2Cell.resize( numberOfInterfaces );

    this->BoundaryConditionList.resize( reader.BCs.size() );
    this->FaceAnalyzerList.resize( reader.FaceAnalyzers.size() );

    // ========================================================================
    //              Read BC data
    // ========================================================================
    for( int currentBC = 0; currentBC < reader.BCs.size(); ++currentBC )
    {
        this->BoundaryConditionList[currentBC] = reader.BCs[currentBC];
    }

    // ========================================================================
    //              Read FA data
    // ========================================================================
    for( int currentFA = 0; currentFA < reader.FaceAnalyzers.size(); ++currentFA )
    {
        this->FaceAnalyzerList[currentFA] = reader.FaceAnalyzers[currentFA];
    }

    // ========================================================================
    //              Read Node data
    // ========================================================================
    for( idType node = 0; node < numberOfNodes; ++node)
    {
        this->NodeCenter[node] = reader.Nodes[node];
    }
    // ========================================================================

    // ========================================================================
    //              Read Cell data
    // ========================================================================
    for( idType cell = 0; cell < numberOfCells; ++cell )
    {
        PrimitiveVariable prim;
        prim.rho = this->fluidParam.referencePrim.rho;
        prim.U   = 0.0;
        prim.V   = 0.0;
        prim.L   = this->fluidParam.referencePrim.L;

        ConservedVariable cons = prim2cons( prim );

        this->CellData   [cell]  = cons;
        this->CellDataOld[cell]  = cons;

        this->Cell2Node[cell] = reader.Cell2Node[cell];

        this->Cell2Interface[cell] = reader.Cell2Face[cell];

        this->CellCenter [cell] = reader.CellCenter [cell];
        this->CellVolume [cell] = reader.CellVolume [cell];
        this->CellMinDx  [cell] = reader.CellMinDx  [cell];
        this->CellLSCoeff[cell] = reader.CellLSCoeff[cell];

        this->CellBoundaryCondition[cell] = reader.Cell2BC[cell];
    }
    // ========================================================================

    // ========================================================================
    //              Read Interface data
    // ========================================================================
    for( idType interface = 0; interface < numberOfInterfaces; ++interface )
    {

        this->Interface2Node[interface][0] = reader.Face2Node[interface][0];
        this->Interface2Node[interface][1] = reader.Face2Node[interface][1];

        this->Interface2Cell[interface][0] = reader.Face2Cell[interface][0];
        this->Interface2Cell[interface][1] = reader.Face2Cell[interface][1];

        this->InterfaceCenter  [interface] = reader.FaceCenter  [interface];
        this->InterfaceNormal  [interface] = reader.FaceNormal  [interface];
        this->InterfaceArea    [interface] = reader.FaceArea    [interface];
        this->InterfaceDistance[interface] = reader.FaceDistance[interface];

        this->InterfaceAdd2Cell[interface][0] = reader.Face2CellAdd[interface][0];
        this->InterfaceAdd2Cell[interface][1] = reader.Face2CellAdd[interface][1];
    }
    // ========================================================================

    int breakPoint = 0;
    return true;
}

void GKSSolverPush::storeDataOld(idType id)
{
    this->CellDataOld[id] = this->CellData[id];
}

void GKSSolverPush::updateCell(const idType id)
{
    // ========================================================================
    //                      Update conservative Variables
    // ========================================================================
    CellData[id].rho  += CellUpdate[id].rho  / CellVolume[id];
    CellData[id].rhoU += CellUpdate[id].rhoU / CellVolume[id];
    CellData[id].rhoV += CellUpdate[id].rhoV / CellVolume[id];
    CellData[id].rhoE += CellUpdate[id].rhoE / CellVolume[id];
    // ========================================================================

    // ========================================================================
    //                      Reset Update storage
    // ========================================================================
    CellUpdate[id].rho  = 0.0;
    CellUpdate[id].rhoU = 0.0;
    CellUpdate[id].rhoV = 0.0;
    CellUpdate[id].rhoE = 0.0;
    // ========================================================================


}

void GKSSolverPush::applyFlux(idType id, ConservedVariable flux)
{
    if( this->InterfaceAdd2Cell[id][0] )
    {
        #pragma omp atomic
        this->CellUpdate[ Interface2Cell[id][0] ].rho  += flux.rho ;
        #pragma omp atomic
        this->CellUpdate[ Interface2Cell[id][0] ].rhoU += flux.rhoU;
        #pragma omp atomic
        this->CellUpdate[ Interface2Cell[id][0] ].rhoV += flux.rhoV;
        #pragma omp atomic
        this->CellUpdate[ Interface2Cell[id][0] ].rhoE += flux.rhoE;
    }
    
    if( this->InterfaceAdd2Cell[id][1] )
    {
        #pragma omp atomic
        this->CellUpdate[ Interface2Cell[id][1] ].rho  -= flux.rho ;
        #pragma omp atomic
        this->CellUpdate[ Interface2Cell[id][1] ].rhoU -= flux.rhoU;
        #pragma omp atomic
        this->CellUpdate[ Interface2Cell[id][1] ].rhoV -= flux.rhoV;
        #pragma omp atomic
        this->CellUpdate[ Interface2Cell[id][1] ].rhoE -= flux.rhoE;
    }

    this->InterfaceFlux[id] = flux;
}

bool GKSSolverPush::isGhostCell(const idType & id)
{
    return this->CellBoundaryCondition[id] != -1;
}

ConservedVariable GKSSolverPush::getCellData(idType id)
{
    return this->CellData[id];
}

ConservedVariable GKSSolverPush::getCellDataOld(idType id)
{
    return this->CellDataOld[id];
}

ConservedVariable GKSSolverPush::getCellGradientX(idType id)
{
    return this->CellGradientX[id];
}

ConservedVariable GKSSolverPush::getCellGradientY(idType id)
{
    return this->CellGradientY[id];
}

idType GKSSolverPush::getCell2Interface(idType cell, idType face)
{
    return this->Cell2Interface[cell][face];
}

Vec2 GKSSolverPush::getCellCenter(idType id)
{
    return this->CellCenter[id];
}

double GKSSolverPush::getCellMinDx(idType id)
{
    return this->CellMinDx[id];
}

array<double, 3> GKSSolverPush::getCellLSCoeff(idType id)
{
    return this->CellLSCoeff[id];
}

int GKSSolverPush::getCellBoundaryCondition(idType id)
{
    return this->CellBoundaryCondition[id];
}

double GKSSolverPush::getInterfaceArea(idType id)
{
    return this->InterfaceArea[id];
}

double GKSSolverPush::getInterfaceDistance(idType id)
{
    return this->InterfaceDistance[id];
}

idType GKSSolverPush::getPosCell(idType id)
{
    return this->Interface2Cell[id][0];
}

idType GKSSolverPush::getNegCell(idType id)
{
    return this->Interface2Cell[id][1];
}

Vec2 GKSSolverPush::getInterfaceNormal(idType id)
{
    return this->InterfaceNormal[id];
}

Vec2 GKSSolverPush::getInterfaceCenter(idType id)
{
    return this->InterfaceCenter[id];
}

void GKSSolverPush::setData(idType id, ConservedVariable cons)
{
    this->CellData[id] = cons;
}

void GKSSolverPush::setCellGradientX(idType id, ConservedVariable dWdx)
{
    this->CellGradientX[id] = dWdx;
}

void GKSSolverPush::setCellGradientY(idType id, ConservedVariable dWdy)
{
    this->CellGradientY[id] = dWdy;
}

Vec2 GKSSolverPush::getNode(idType node)
{
    return this->NodeCenter[node];
}

idType GKSSolverPush::getNumberOfNodes()
{
    return this->numberOfNodes;
}

idType GKSSolverPush::getNumberOfCells()
{
    return this->numberOfCells;
}

idType GKSSolverPush::getNumberOfInterfaces()
{
    return this->numberOfInterfaces;
}

array<idType, 4> GKSSolverPush::getCell2Node(idType cell)
{
    return this->Cell2Node[cell];
}


