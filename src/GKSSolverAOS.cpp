
#include "GKSSolverAOS.h"
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

GKSSolverAOS::GKSSolverAOS() : GKSSolver()
{
}

GKSSolverAOS::GKSSolverAOS(Parameters param, FluidParameter fluidParam)
         : GKSSolver(param, fluidParam)
{
}

GKSSolverAOS::~GKSSolverAOS()
{
}

void GKSSolverAOS::readMeshFromMeshObject(const GKSMesh& origin)
{
    this->numberOfNodes        = origin.NodeList.size();

    this->numberOfCells        = origin.CellList.size();

    this->numberOfInterfaces   = origin.InterfaceList.size();

    // ========================================================================
    //              Allocate Memory
    // ========================================================================

    this->Cells.resize( numberOfCells );

    this->Interfaces.resize( numberOfInterfaces );

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

    // ========================================================================
    //              Read Cell data
    // ========================================================================
    for( Cell* currentCell : origin.CellList )
    {
        this->Cells[currentCell->getID()-1].Data    = currentCell->getCons();
        this->Cells[currentCell->getID()-1].DataOld = currentCell->getConsOld();

        for( int i = 0; i < 4; i++)
        {
            if(currentCell->getInterface(i) != NULL)
                this->Cells[currentCell->getID()-1].Interfaces[i] = currentCell->getInterface(i)->getID()-1;
            else
                this->Cells[currentCell->getID()-1].Interfaces[i] = -1;
        }

        this->Cells[currentCell->getID()-1].Center = currentCell->getCenter();
        this->Cells[currentCell->getID()-1].Volume = currentCell->getVolume();
        this->Cells[currentCell->getID()-1].MinDx  = currentCell->getMinDx();
    }
    // ========================================================================

    // ========================================================================
    //              Read Interface data
    // ========================================================================
    for( Interface* currentInterface : origin.InterfaceList )
    {
        this->Interfaces[currentInterface->getID()-1].Flux = currentInterface->getTimeIntegratedFlux();

        for( int i = 0; i < 2; i++)
            this->Interfaces[currentInterface->getID()-1].Cells[i] = currentInterface->getCell(i)->getID()-1;

        this->Interfaces[currentInterface->getID()-1].Normal   = currentInterface->getNormal();
        this->Interfaces[currentInterface->getID()-1].Area     = currentInterface->getArea();
        this->Interfaces[currentInterface->getID()-1].Distance = currentInterface->getDistance2CellCenter(0) + currentInterface->getDistance2CellCenter(1);
    }
    // ========================================================================

    int breakPoint = 0;
}

void GKSSolverAOS::writeDataToMeshObject(const GKSMesh & target)
{
    for( Cell* currentCell : target.CellList )
        currentCell->setCons( this->Cells[ currentCell->getID()-1 ].Data );
}

void GKSSolverAOS::storeDataOld(idType id)
{
    this->Cells[id].DataOld = this->Cells[id].Data;
}

void GKSSolverAOS::updateCell(const idType id)
{
    // ========================================================================
    //                      Compute Flux signs
    // ========================================================================
    double fluxSign[4];
    fluxSign[0] = ( Interfaces[ Cells[id].Interfaces[0] ].Cells[0] == id ) ? (1.0) : (-1.0);
    fluxSign[1] = ( Interfaces[ Cells[id].Interfaces[1] ].Cells[0] == id ) ? (1.0) : (-1.0);
    fluxSign[2] = ( Interfaces[ Cells[id].Interfaces[2] ].Cells[0] == id ) ? (1.0) : (-1.0);
    fluxSign[3] = ( Interfaces[ Cells[id].Interfaces[3] ].Cells[0] == id ) ? (1.0) : (-1.0);
    // ========================================================================

    // ========================================================================
    //                      Update conservative Variables
    // ========================================================================
    Cells[id].Data.rho  += ( fluxSign[0] * Interfaces[ Cells[id].Interfaces[0] ].Flux.rho
                           + fluxSign[1] * Interfaces[ Cells[id].Interfaces[1] ].Flux.rho
                           + fluxSign[2] * Interfaces[ Cells[id].Interfaces[2] ].Flux.rho
                           + fluxSign[3] * Interfaces[ Cells[id].Interfaces[3] ].Flux.rho
                           ) / Cells[id].Volume;

    Cells[id].Data.rhoU += ( fluxSign[0] * Interfaces[ Cells[id].Interfaces[0] ].Flux.rhoU
                           + fluxSign[1] * Interfaces[ Cells[id].Interfaces[1] ].Flux.rhoU
                           + fluxSign[2] * Interfaces[ Cells[id].Interfaces[2] ].Flux.rhoU
                           + fluxSign[3] * Interfaces[ Cells[id].Interfaces[3] ].Flux.rhoU
                           ) / Cells[id].Volume;

    Cells[id].Data.rhoV += ( fluxSign[0] * Interfaces[ Cells[id].Interfaces[0] ].Flux.rhoV
                           + fluxSign[1] * Interfaces[ Cells[id].Interfaces[1] ].Flux.rhoV
                           + fluxSign[2] * Interfaces[ Cells[id].Interfaces[2] ].Flux.rhoV
                           + fluxSign[3] * Interfaces[ Cells[id].Interfaces[3] ].Flux.rhoV
                           ) / Cells[id].Volume;

    Cells[id].Data.rhoE += ( fluxSign[0] * Interfaces[ Cells[id].Interfaces[0] ].Flux.rhoE
                           + fluxSign[1] * Interfaces[ Cells[id].Interfaces[1] ].Flux.rhoE
                           + fluxSign[2] * Interfaces[ Cells[id].Interfaces[2] ].Flux.rhoE
                           + fluxSign[3] * Interfaces[ Cells[id].Interfaces[3] ].Flux.rhoE
                           ) / Cells[id].Volume;
    // ========================================================================
}

void GKSSolverAOS::applyFlux(idType id, ConservedVariable flux)
{
    this->Interfaces[id].Flux = flux;
}

bool GKSSolverAOS::isGhostCell(const idType & id)
{
    return  this->Cells[id].Interfaces[0] == -1
         || this->Cells[id].Interfaces[1] == -1
         || this->Cells[id].Interfaces[2] == -1
         || this->Cells[id].Interfaces[3] == -1;
}

ConservedVariable GKSSolverAOS::getCellData(idType id)
{
    return this->Cells[id].Data;
}

ConservedVariable GKSSolverAOS::getCellDataOld(idType id)
{
    return this->Cells[id].DataOld;
}

double GKSSolverAOS::getCellMinDx(idType id)
{
    return this->Cells[id].MinDx;
}

double GKSSolverAOS::getInterfaceArea(idType id)
{
    return this->Interfaces[id].Area;
}

double GKSSolverAOS::getInterfaceDistance(idType id)
{
    return this->Interfaces[id].Distance;
}

idType GKSSolverAOS::getPosCell(idType id)
{
    return this->Interfaces[id].Cells[0];
}

idType GKSSolverAOS::getNegCell(idType id)
{
    return this->Interfaces[id].Cells[1];
}

Vec2 GKSSolverAOS::getInterfaceNormal(idType id)
{
    return this->Interfaces[id].Normal;
}

void GKSSolverAOS::setData(idType id, ConservedVariable cons)
{
    this->Cells[id].Data = cons;
}


