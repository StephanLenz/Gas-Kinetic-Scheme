// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      GKSSolverPush.h
//
//      Function:
//          Generation and Storage of mesh
//          Control of simulation
//          Data Analysis
//          File Output
//
// ============================================================================

#ifndef GKSSolverPush_H
#define GKSSolverPush_H

#include "GKSSolver.h"
#include "BoundaryCondition.h"
#include "Types.h"
#include <vector>
#include <array>
#include <list>
#include <string>
#include <chrono>

using namespace std;

class GKSSolverPush : public GKSSolver
{
// ====================================================================================================================
// ====================================================================================================================
//                      Attributes
// ====================================================================================================================
// ====================================================================================================================
private:

    // ========================================================================
    //              data
    // ========================================================================

    vector<ConservedVariable> CellData;
    vector<ConservedVariable> CellDataOld;

    vector<ConservedVariable> CellUpdate;

    vector<ConservedVariable> CellGradientX;
    vector<ConservedVariable> CellGradientY;

    vector<ConservedVariable> InterfaceFlux;

    // ========================================================================
    //              Connectivity
    // ========================================================================

    vector< array<idType, 4> > Cell2Node;
    vector< array<idType, 4> > Cell2Interface;
    vector< idType >           CellBoundaryCondition;

    vector< array<idType, 2> > Interface2Node;
    vector< array<idType, 2> > Interface2Cell;

    vector< array<bool, 2> > InterfaceAdd2Cell;

    // ========================================================================
    //              Geometry
    // ========================================================================

    vector<Vec2> NodeCenter;

    vector<Vec2>              CellCenter;
    vector<double>            CellVolume;
    vector<double>            CellMinDx;
    vector< array<double,3> > CellLSCoeff;

    vector<Vec2> InterfaceCenter;
    vector<Vec2> InterfaceNormal;
    vector<double> InterfaceDistance;
    vector<double> InterfaceArea;
    vector< array<double,2> > Interface2CellCenterDistance;

// ====================================================================================================================
// ====================================================================================================================
//                      Methods
// ====================================================================================================================
// ====================================================================================================================
public:
	GKSSolverPush();

    GKSSolverPush(Parameters param, FluidParameter fluidParam);

	~GKSSolverPush();

    // ========================================================================
    //              Communication methods
    // ========================================================================

    virtual bool readProblem( string filename );

    // ========================================================================
    //              Simulation Control
    // ========================================================================

    virtual void storeDataOld(idType id);

    virtual void updateCell(const idType id);

    // ========================================================================
    //              Flux computation subroutines
    // ========================================================================

    virtual void applyFlux(idType id, ConservedVariable flux);

    // ========================================================================
    //              Util
    // ========================================================================

    virtual bool isGhostCell(const idType& id);

    virtual ConservedVariable getCellData(idType id);
    virtual ConservedVariable getCellDataOld(idType id);

    virtual ConservedVariable getCellGradientX(idType id);
    virtual ConservedVariable getCellGradientY(idType id);

    virtual idType getCell2Interface( idType cell, idType face );
    
    virtual Vec2   getCellCenter(idType id);
    virtual double getCellMinDx(idType id);
    virtual array<double,3> getCellLSCoeff(idType id);

    virtual int getCellBoundaryCondition(idType id);

    virtual double getInterfaceArea(idType id);
    virtual double getInterfaceDistance(idType id);

    virtual idType getPosCell(idType id);
    virtual idType getNegCell(idType id);

    virtual Vec2 getInterfaceNormal(idType id);
    virtual Vec2 getInterfaceCenter(idType id);

    virtual void setData(idType id, ConservedVariable cons);

    virtual void setCellGradientX(idType id, ConservedVariable dWdx);
    virtual void setCellGradientY(idType id, ConservedVariable dWdy);

    virtual Vec2 getNode(idType node);

    virtual array<idType,4> getCell2Node(idType cell);
};

#endif