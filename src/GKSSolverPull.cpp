
#include "GKSSolverPull.h"
#include "Cell.h"
#include "Interface.h"
#include "BoundaryCondition.h"
#include "InterfaceBC.h"
#include "Types.h"
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

GKSSolverPull::GKSSolverPull() : GKSSolver()
{
}

GKSSolverPull::GKSSolverPull(Parameters param, FluidParameter fluidParam)
         : GKSSolver(param, fluidParam)
{
}

GKSSolverPull::~GKSSolverPull()
{
}

void GKSSolverPull::readMeshFromMeshObject(const GKSMesh& origin)
{
    this->numberOfNodes        = origin.NodeList.size();

    this->numberOfCells        = origin.CellList.size();

    this->numberOfInterfaces   = origin.InterfaceList.size();

    // ========================================================================
    //              Allocate Memory
    // ========================================================================

    this->NodeCenter.resize( numberOfNodes );

    this->CellData.resize( numberOfCells );
    this->CellDataOld.resize( numberOfCells );

    this->Cell2Node.resize( numberOfCells );
    this->Cell2Interface.resize( numberOfCells );
    this->CellBoundaryCondition.resize( numberOfCells );

    this->CellCenter.resize( numberOfCells );
    this->CellVolume.resize( numberOfCells );
    this->CellMinDx.resize( numberOfCells );

    this->InterfaceFlux.resize( numberOfInterfaces );

    this->Interface2Node.resize( numberOfInterfaces );
    this->Interface2Cell.resize( numberOfInterfaces );

    this->InterfaceCenter.resize( numberOfInterfaces );
    this->InterfaceNormal.resize( numberOfInterfaces );
    this->InterfaceArea.resize( numberOfInterfaces );
    this->InterfaceDistance.resize( numberOfInterfaces );
    this->Interface2CellCenterDistance.resize( numberOfInterfaces );


    // ========================================================================
    //              Read BC data
    // ========================================================================
    for( BoundaryCondition* currentBC: origin.BoundaryConditionList )
    {
        this->BoundaryConditionList.push_back(*currentBC);
        for(Cell* currentCell : origin.CellList)
        {
            if( currentBC == currentCell->getBoundaryConditionPointer() )
            {
                this->BoundaryConditionList.back().addCell( currentCell->getID() - 1 );
                this->BoundaryConditionList.back().addNeighborCell( currentCell->findNeighborInDomain()->getID() - 1 );
            }
        }
    }

    // ========================================================================
    //              Read Node data
    // ========================================================================
    for( Node* currentNode : origin.NodeList )
    {
        this->NodeCenter[currentNode->ID-1] = *currentNode;
    }
    // ========================================================================

    // ========================================================================
    //              Read Cell data
    // ========================================================================
    for( Cell* currentCell : origin.CellList )
    {
        this->CellData   [currentCell->getID()-1] = currentCell->getCons();
        this->CellDataOld[currentCell->getID()-1] = currentCell->getConsOld();

        for( int i = 0; i < 4; i++)
            this->Cell2Node[currentCell->getID()-1][i] = currentCell->getNode(i).ID-1;

        for( int i = 0; i < 4; i++)
        {
            if(currentCell->getInterface(i) != NULL)
                this->Cell2Interface[currentCell->getID()-1][i] = currentCell->getInterface(i)->getID()-1;
            else
                this->Cell2Interface[currentCell->getID()-1][i] = -1;
        }

        this->CellCenter[currentCell->getID()-1] = currentCell->getCenter();
        this->CellVolume[currentCell->getID()-1] = currentCell->getVolume();
        this->CellMinDx [currentCell->getID()-1] = currentCell->getMinDx();

        this->CellBoundaryCondition[currentCell->getID()-1] = -1;
        for( int i = 0; i < origin.BoundaryConditionList.size(); ++i)
        {
            if( origin.BoundaryConditionList[i] == currentCell->getBoundaryConditionPointer() )
            {
                this->CellBoundaryCondition[currentCell->getID()-1] = i;
                break;
            }
        }
    }
    // ========================================================================

    // ========================================================================
    //              Read Interface data
    // ========================================================================
    for( Interface* currentInterface : origin.InterfaceList )
    {
        this->InterfaceFlux[currentInterface->getID()-1] = currentInterface->getTimeIntegratedFlux();

        for( int i = 0; i < 2; i++)
            this->Interface2Node[currentInterface->getID()-1][i] = currentInterface->getNode(i)->ID-1;

        for( int i = 0; i < 2; i++)
            this->Interface2Cell[currentInterface->getID()-1][i] = currentInterface->getCell(i)->getID()-1;

        this->InterfaceCenter[currentInterface->getID()-1]   = currentInterface->getCenter();
        this->InterfaceNormal[currentInterface->getID()-1]   = currentInterface->getNormal();
        this->InterfaceArea[currentInterface->getID()-1]     = currentInterface->getArea();
        this->InterfaceDistance[currentInterface->getID()-1] = currentInterface->getDistance2CellCenter(0) + currentInterface->getDistance2CellCenter(1);

        for( int i = 0; i < 2; i++)
            this->Interface2CellCenterDistance[currentInterface->getID()-1][i] = currentInterface->getDistance2CellCenter(i);
    }
    // ========================================================================

    int breakPoint = 0;
}

void GKSSolverPull::writeDataToMeshObject(const GKSMesh & target)
{
    for( Cell* currentCell : target.CellList )
        currentCell->setCons( this->CellData[ currentCell->getID()-1 ] );
}

void GKSSolverPull::storeDataOld(idType id)
{
    this->CellDataOld[id] = this->CellData[id];
}

void GKSSolverPull::updateCell(const idType id)
{
    // ========================================================================
    //                      Compute Flux signs
    // ========================================================================
    double fluxSign[4];
    fluxSign[0] = ( Interface2Cell[ Cell2Interface[id][0] ][0] == id ) ? (1.0) : (-1.0);
    fluxSign[1] = ( Interface2Cell[ Cell2Interface[id][1] ][0] == id ) ? (1.0) : (-1.0);
    fluxSign[2] = ( Interface2Cell[ Cell2Interface[id][2] ][0] == id ) ? (1.0) : (-1.0);
    fluxSign[3] = ( Interface2Cell[ Cell2Interface[id][3] ][0] == id ) ? (1.0) : (-1.0);
    // ========================================================================

    // ========================================================================
    //                      Update conservative Variables
    // ========================================================================
    CellData[id].rho  += ( fluxSign[0] * InterfaceFlux[ Cell2Interface[id][0] ].rho
                         + fluxSign[1] * InterfaceFlux[ Cell2Interface[id][1] ].rho
                         + fluxSign[2] * InterfaceFlux[ Cell2Interface[id][2] ].rho
                         + fluxSign[3] * InterfaceFlux[ Cell2Interface[id][3] ].rho
                         ) / CellVolume[id];

    CellData[id].rhoU += ( fluxSign[0] * InterfaceFlux[ Cell2Interface[id][0] ].rhoU
                         + fluxSign[1] * InterfaceFlux[ Cell2Interface[id][1] ].rhoU
                         + fluxSign[2] * InterfaceFlux[ Cell2Interface[id][2] ].rhoU
                         + fluxSign[3] * InterfaceFlux[ Cell2Interface[id][3] ].rhoU
                         ) / CellVolume[id];

    CellData[id].rhoV += ( fluxSign[0] * InterfaceFlux[ Cell2Interface[id][0] ].rhoV
                         + fluxSign[1] * InterfaceFlux[ Cell2Interface[id][1] ].rhoV
                         + fluxSign[2] * InterfaceFlux[ Cell2Interface[id][2] ].rhoV
                         + fluxSign[3] * InterfaceFlux[ Cell2Interface[id][3] ].rhoV
                         ) / CellVolume[id];

    CellData[id].rhoE += ( fluxSign[0] * InterfaceFlux[ Cell2Interface[id][0] ].rhoE
                         + fluxSign[1] * InterfaceFlux[ Cell2Interface[id][1] ].rhoE
                         + fluxSign[2] * InterfaceFlux[ Cell2Interface[id][2] ].rhoE
                         + fluxSign[3] * InterfaceFlux[ Cell2Interface[id][3] ].rhoE
                         ) / CellVolume[id];
    // ========================================================================
}

void GKSSolverPull::applyFlux(idType id, ConservedVariable flux)
{
    this->InterfaceFlux[id] = flux;
}

bool GKSSolverPull::isGhostCell(const idType & id)
{
    return this->CellBoundaryCondition[id] != -1;
}

ConservedVariable GKSSolverPull::getCellData(idType id)
{
    return this->CellData[id];
}

ConservedVariable GKSSolverPull::getCellDataOld(idType id)
{
    return this->CellDataOld[id];
}

double GKSSolverPull::getCellMinDx(idType id)
{
    return this->CellMinDx[id];
}

double GKSSolverPull::getInterfaceArea(idType id)
{
    return this->InterfaceArea[id];
}

double GKSSolverPull::getInterfaceDistance(idType id)
{
    return this->InterfaceDistance[id];
}

idType GKSSolverPull::getPosCell(idType id)
{
    return this->Interface2Cell[id][0];
}

idType GKSSolverPull::getNegCell(idType id)
{
    return this->Interface2Cell[id][1];
}

Vec2 GKSSolverPull::getInterfaceNormal(idType id)
{
    return this->InterfaceNormal[id];
}

void GKSSolverPull::setData(idType id, ConservedVariable cons)
{
    this->CellData[id] = cons;
}


