// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      GKSSolverSOA.h
//
//      Function:
//          Generation and Storage of mesh
//          Control of simulation
//          Data Analysis
//          File Output
//
// ============================================================================

#ifndef GKSSolverSOA_H
#define GKSSolverSOA_H

#include "GKSSolver.h"
#include "GKSMesh.h"
#include "Cell.h"
#include "Interface.h"
#include "BoundaryCondition.h"
#include "InterfaceBC.h"
#include "Types.h"
#include <vector>
#include <array>
#include <list>
#include <string>
#include <chrono>

using namespace std;

class GKSSolverSOA : public GKSSolver
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

    vector<double> rho;
    vector<double> rhoU;
    vector<double> rhoV;
    vector<double> rhoE;

    vector<double> rho_Old;
    vector<double> rhoU_Old;
    vector<double> rhoV_Old;
    vector<double> rhoE_Old;
    
    vector<double> F_rho;
    vector<double> F_rhoU;
    vector<double> F_rhoV;
    vector<double> F_rhoE;

    // ========================================================================
    //              Connectivity
    // ========================================================================
    
    vector<idType> Cell2Interface_0;
    vector<idType> Cell2Interface_1;
    vector<idType> Cell2Interface_2;
    vector<idType> Cell2Interface_3;

    vector<idType> Interface2CellPos;
    vector<idType> Interface2CellNeg;

    // ========================================================================
    //              Geometry
    // ========================================================================

    vector<double> CellVolume;
    vector<double> CellMinDx;

    vector<double> InterfaceNormalX;
    vector<double> InterfaceNormalY;

    vector<double> InterfaceDistance;
    vector<double> InterfaceArea;

// ====================================================================================================================
// ====================================================================================================================
//                      Methods
// ====================================================================================================================
// ====================================================================================================================
public:
	GKSSolverSOA();

    GKSSolverSOA(Parameters param, FluidParameter fluidParam);

	~GKSSolverSOA();

    // ========================================================================
    //              Communication methods
    // ========================================================================

    virtual void readMeshFromMeshObject( const GKSMesh& origin );

    virtual void writeDataToMeshObject( const GKSMesh& target );

    // ========================================================================
    //              Simulation Control
    // ========================================================================

    virtual void updateCell(const idType id);

    // ========================================================================
    //              Flux computation subroutines
    // ========================================================================

    virtual void storeDataOld(idType id);

    virtual void applyFlux(idType id, ConservedVariable flux);

    // ========================================================================
    //              Util
    // ========================================================================

    virtual bool isGhostCell(const idType& id);

    virtual ConservedVariable getCellData(idType id);
    virtual ConservedVariable getCellDataOld(idType id);

    virtual double getCellMinDx(idType id);

    virtual double getInterfaceArea(idType id);
    virtual double getInterfaceDistance(idType id);

    virtual idType getPosCell(idType id);
    virtual idType getNegCell(idType id);

    virtual Vec2 getInterfaceNormal(idType id);

    virtual void setData(idType id, ConservedVariable cons);
};

#endif