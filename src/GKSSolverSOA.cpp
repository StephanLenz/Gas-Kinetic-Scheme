
#include "GKSSolverSOA.h"
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

GKSSolverSOA::GKSSolverSOA() : GKSSolver()
{
}

GKSSolverSOA::GKSSolverSOA(Parameters param, FluidParameter fluidParam)
         : GKSSolver(param, fluidParam)
{
}

GKSSolverSOA::~GKSSolverSOA()
{
}

void GKSSolverSOA::readMeshFromMeshObject(const GKSMesh& origin)
{
    this->numberOfNodes        = origin.NodeList.size();

    this->numberOfCells        = origin.CellList.size();

    this->numberOfInterfaces   = origin.InterfaceList.size();

    // ========================================================================
    //              Allocate Memory
    // ========================================================================

    rho.resize( numberOfCells );
    rhoU.resize( numberOfCells );
    rhoV.resize( numberOfCells );
    rhoE.resize( numberOfCells );

    rho_Old.resize( numberOfCells );
    rhoU_Old.resize( numberOfCells );
    rhoV_Old.resize( numberOfCells );
    rhoE_Old.resize( numberOfCells );
    
    Cell2Interface_0.resize( numberOfCells );
    Cell2Interface_1.resize( numberOfCells );
    Cell2Interface_2.resize( numberOfCells );
    Cell2Interface_3.resize( numberOfCells );

    CellVolume.resize( numberOfCells );
    CellMinDx.resize( numberOfCells );
    
    F_rho.resize( numberOfInterfaces );
    F_rhoU.resize( numberOfInterfaces );
    F_rhoV.resize( numberOfInterfaces );
    F_rhoE.resize( numberOfInterfaces );
    
    InterfaceNormalX.resize( numberOfInterfaces );
    InterfaceNormalY.resize( numberOfInterfaces );

    Interface2CellPos.resize( numberOfInterfaces );
    Interface2CellNeg.resize( numberOfInterfaces );

    InterfaceDistance.resize( numberOfInterfaces );
    InterfaceArea.resize( numberOfInterfaces );

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
        this->rho [currentCell->getID()-1] = currentCell->getCons().rho;
        this->rhoU[currentCell->getID()-1] = currentCell->getCons().rhoU;
        this->rhoV[currentCell->getID()-1] = currentCell->getCons().rhoV;
        this->rhoE[currentCell->getID()-1] = currentCell->getCons().rhoE;

        this->rho_Old [currentCell->getID()-1] = currentCell->getConsOld().rho;
        this->rhoU_Old[currentCell->getID()-1] = currentCell->getConsOld().rhoU;
        this->rhoV_Old[currentCell->getID()-1] = currentCell->getConsOld().rhoV;
        this->rhoE_Old[currentCell->getID()-1] = currentCell->getConsOld().rhoE;

        if(currentCell->getInterface(0) != NULL)
            this->Cell2Interface_0[currentCell->getID()-1] = currentCell->getInterface(0)->getID()-1;
        else
            this->Cell2Interface_0[currentCell->getID()-1] = -1;
        
        if(currentCell->getInterface(1) != NULL)
            this->Cell2Interface_1[currentCell->getID()-1] = currentCell->getInterface(1)->getID()-1;
        else
            this->Cell2Interface_1[currentCell->getID()-1] = -1;
        
        if(currentCell->getInterface(2) != NULL)
            this->Cell2Interface_2[currentCell->getID()-1] = currentCell->getInterface(2)->getID()-1;
        else
            this->Cell2Interface_2[currentCell->getID()-1] = -1;
        
        if(currentCell->getInterface(3) != NULL)
            this->Cell2Interface_3[currentCell->getID()-1] = currentCell->getInterface(3)->getID()-1;
        else
            this->Cell2Interface_3[currentCell->getID()-1] = -1;

        this->CellVolume[currentCell->getID()-1] = currentCell->getVolume();
        this->CellMinDx [currentCell->getID()-1] = currentCell->getMinDx();
    }
    // ========================================================================

    // ========================================================================
    //              Read Interface data
    // ========================================================================
    for( Interface* currentInterface : origin.InterfaceList )
    {
        this->F_rho [currentInterface->getID()-1] = currentInterface->getTimeIntegratedFlux().rho;
        this->F_rhoU[currentInterface->getID()-1] = currentInterface->getTimeIntegratedFlux().rhoU;
        this->F_rhoV[currentInterface->getID()-1] = currentInterface->getTimeIntegratedFlux().rhoV;
        this->F_rhoE[currentInterface->getID()-1] = currentInterface->getTimeIntegratedFlux().rhoE;

        this->Interface2CellPos[currentInterface->getID()-1] = currentInterface->getCell(0)->getID()-1;
        this->Interface2CellNeg[currentInterface->getID()-1] = currentInterface->getCell(1)->getID()-1;

        this->InterfaceNormalX[currentInterface->getID()-1] = currentInterface->getNormal().x;
        this->InterfaceNormalY[currentInterface->getID()-1] = currentInterface->getNormal().y;

        this->InterfaceArea[currentInterface->getID()-1]     = currentInterface->getArea();
        this->InterfaceDistance[currentInterface->getID()-1] = currentInterface->getDistance2CellCenter(0) + currentInterface->getDistance2CellCenter(1);
    }
    // ========================================================================
}

void GKSSolverSOA::writeDataToMeshObject(const GKSMesh & target)
{
    for( Cell* currentCell : target.CellList )
    {
        ConservedVariable tmp;
        tmp.rho  = rho [ currentCell->getID()-1 ];
        tmp.rhoU = rhoU[ currentCell->getID()-1 ];
        tmp.rhoV = rhoV[ currentCell->getID()-1 ];
        tmp.rhoE = rhoE[ currentCell->getID()-1 ];
        currentCell->setCons( tmp );
    }
}

void GKSSolverSOA::updateCell(const idType id)
{
    // ========================================================================
    //                      Compute Flux signs
    // ========================================================================
    double fluxSign[4];
    fluxSign[0] = ( Interface2CellPos[ Cell2Interface_0[id] ] == id ) ? (1.0) : (-1.0);
    fluxSign[1] = ( Interface2CellPos[ Cell2Interface_1[id] ] == id ) ? (1.0) : (-1.0);
    fluxSign[2] = ( Interface2CellPos[ Cell2Interface_2[id] ] == id ) ? (1.0) : (-1.0);
    fluxSign[3] = ( Interface2CellPos[ Cell2Interface_3[id] ] == id ) ? (1.0) : (-1.0);
    // ========================================================================

    // ========================================================================
    //                      Update conservative Variables
    // ========================================================================
    rho [id] += ( fluxSign[0] * F_rho [ Cell2Interface_0[id] ]
                + fluxSign[1] * F_rho [ Cell2Interface_1[id] ]
                + fluxSign[2] * F_rho [ Cell2Interface_2[id] ]
                + fluxSign[3] * F_rho [ Cell2Interface_3[id] ]
                ) / CellVolume[id];

    rhoU[id] += ( fluxSign[0] * F_rhoU[ Cell2Interface_0[id] ]
                + fluxSign[1] * F_rhoU[ Cell2Interface_1[id] ]
                + fluxSign[2] * F_rhoU[ Cell2Interface_2[id] ]
                + fluxSign[3] * F_rhoU[ Cell2Interface_3[id] ]
                ) / CellVolume[id];

    rhoV[id] += ( fluxSign[0] * F_rhoV[ Cell2Interface_0[id] ]
                + fluxSign[1] * F_rhoV[ Cell2Interface_1[id] ]
                + fluxSign[2] * F_rhoV[ Cell2Interface_2[id] ]
                + fluxSign[3] * F_rhoV[ Cell2Interface_3[id] ]
                ) / CellVolume[id];

    rhoE[id] += ( fluxSign[0] * F_rhoE[ Cell2Interface_0[id] ]
                + fluxSign[1] * F_rhoE[ Cell2Interface_1[id] ]
                + fluxSign[2] * F_rhoE[ Cell2Interface_2[id] ]
                + fluxSign[3] * F_rhoE[ Cell2Interface_3[id] ]
                ) / CellVolume[id];
    // ========================================================================

    int breakPoint = 0;
}

void GKSSolverSOA::storeDataOld(idType id)
{
    this->rho_Old [id] = this->rho [id];
    this->rhoU_Old[id] = this->rhoU[id];
    this->rhoV_Old[id] = this->rhoV[id];
    this->rhoE_Old[id] = this->rhoE[id];
}

void GKSSolverSOA::applyFlux(idType id, ConservedVariable flux)
{
    this->F_rho [id] = flux.rho ;
    this->F_rhoU[id] = flux.rhoU;
    this->F_rhoV[id] = flux.rhoV;
    this->F_rhoE[id] = flux.rhoE;
}

bool GKSSolverSOA::isGhostCell(const idType & id)
{
    return  this->Cell2Interface_0[id] == -1
         || this->Cell2Interface_1[id] == -1
         || this->Cell2Interface_2[id] == -1
         || this->Cell2Interface_3[id] == -1;
}

ConservedVariable GKSSolverSOA::getCellData(idType id)
{
    ConservedVariable tmp;
    tmp.rho  = rho [id];
    tmp.rhoU = rhoU[id];
    tmp.rhoV = rhoV[id];
    tmp.rhoE = rhoE[id];
    return tmp;
}

ConservedVariable GKSSolverSOA::getCellDataOld(idType id)
{
    ConservedVariable tmp;
    tmp.rho  = rho_Old [id];
    tmp.rhoU = rhoU_Old[id];
    tmp.rhoV = rhoV_Old[id];
    tmp.rhoE = rhoE_Old[id];
    return tmp;
}

double GKSSolverSOA::getCellMinDx(idType id)
{
    return this->CellMinDx[id];
}

double GKSSolverSOA::getInterfaceArea(idType id)
{
    return this->InterfaceArea[id];
}

double GKSSolverSOA::getInterfaceDistance(idType id)
{
    return this->InterfaceDistance[id];
}

idType GKSSolverSOA::getPosCell(idType id)
{
    return this->Interface2CellPos[id];
}

idType GKSSolverSOA::getNegCell(idType id)
{
    return this->Interface2CellNeg[id];
}

Vec2 GKSSolverSOA::getInterfaceNormal(idType id)
{
    Vec2 tmp;
    tmp.x = this->InterfaceNormalX[id];
    tmp.y = this->InterfaceNormalY[id];
    return tmp;
}

void GKSSolverSOA::setData(idType id, ConservedVariable cons)
{
    this->rho [id] = cons.rho ;
    this->rhoU[id] = cons.rhoU;
    this->rhoV[id] = cons.rhoV;
    this->rhoE[id] = cons.rhoE;
}


