// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      GKSSolverAOS.h
//
//      Function:
//          Generation and Storage of mesh
//          Control of simulation
//          Data Analysis
//          File Output
//
// ============================================================================

#ifndef GKSSolverAOS_H
#define GKSSolverAOS_H

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

class GKSSolverAOS : public GKSSolver
{
// ====================================================================================================================
// ====================================================================================================================
//                      Structure Definitions
// ====================================================================================================================
// ====================================================================================================================
private:
    struct CellStruct 
    {
        ConservedVariable Data;         // 4 double = 32 Byte
        ConservedVariable DataOld;      // 4 double = 32 Byte

        Vec2              Center;       // 2 double = 16 Byte
        double            Volume;       // 1 double =  8 Byte
        double            MinDx;        // 1 double =  8 Byte

        array<idType, 4>  Interfaces;   // 4 idType = 16 Byte
                                        //            -------
                                        //            64 + 48 Byte
    };
    // ========================================================================
    struct InterfaceStruct
    {
        ConservedVariable Flux;         // 4 double = 32 Byte
        Vec2              Normal;       // 2 double = 16 Byte
        double            Distance;     // 1 double =  8 Byte
        double            Area;         // 1 double =  8 Byte
        array<idType, 2>  Cells;        // 2 idType =  8 Byte
                                        //            -------
                                        //            64 +  8 Byte
    };
// ====================================================================================================================
// ====================================================================================================================
//                      Attributes
// ====================================================================================================================
// ====================================================================================================================
private:

    vector<CellStruct>      Cells;
    vector<InterfaceStruct> Interfaces;

// ====================================================================================================================
// ====================================================================================================================
//                      Methods
// ====================================================================================================================
// ====================================================================================================================
public:
	GKSSolverAOS();

    GKSSolverAOS(Parameters param, FluidParameter fluidParam);

	~GKSSolverAOS();

    // ========================================================================
    //              Communication methods
    // ========================================================================

    virtual void readMeshFromMeshObject( const GKSMesh& origin );

    virtual void writeDataToMeshObject( const GKSMesh& target );

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

    virtual double getCellMinDx(idType id);

    virtual double getInterfaceArea(idType id);
    virtual double getInterfaceDistance(idType id);

    virtual idType getPosCell(idType id);
    virtual idType getNegCell(idType id);

    virtual Vec2 getInterfaceNormal(idType id);

    virtual void setData(idType id, ConservedVariable cons);
};

#endif