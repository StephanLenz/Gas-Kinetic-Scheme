// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      GKSSolver.h
//
//      Function:
//          Generation and Storage of mesh
//          Control of simulation
//          Data Analysis
//          File Output
//
// ============================================================================

#ifndef GKSSOLVER_H
#define GKSSOLVER_H

#include "GKSMesh.h"
#include "Cell.h"
#include "Interface.h"
#include "BoundaryCondition.h"
#include "InterfaceBC.h"
#include "Types.h"
#include <vector>
#include <list>
#include <string>
#include <chrono>

using namespace std;

class GKSSolver
{
// ====================================================================================================================
// ====================================================================================================================
//                      Attributes
// ====================================================================================================================
// ====================================================================================================================
private:
    idType numberOfNodes;
    idType numberOfCells;
    idType numberOfInterfaces;

    // ========================================================================
    //              data
    // ========================================================================

    vector<ConservedVariable> CellData;
    vector<ConservedVariable> CellDataOld;

    vector<ConservedVariable> InterfaceFlux;

    vector<BoundaryCondition> BoundaryConditionList;

    // ========================================================================
    //              Connectivity
    // ========================================================================

    vector< array<idType, 4> > Cell2Node;
    vector< array<idType, 4> > Cell2Interface;
    vector< idType >           CellBoundaryCondition;

    vector< array<idType, 2> > Interface2Node;
    vector< array<idType, 2> > Interface2Cell;

    // ========================================================================
    //              Geometry
    // ========================================================================

    vector<Vec2> NodeCenter;

    vector<Vec2> CellCenter;
    vector<double> CellVolume;
    vector<double> CellMinDx;

    vector<Vec2> InterfaceCenter;
    vector<Vec2> InterfaceNormal;
    vector<double> InterfaceDistance;
    vector<double> InterfaceArea;
    vector< array<double,2> > Interface2CellCenterDistance;

    // ========================================================================
    //              Fluid Parameters
    // ========================================================================

    FluidParameter fluidParam;
    
    // ========================================================================
    //              Simulation Parameters/Data
    // ========================================================================

    Parameters param;
    double dt;
    double time;

    unsigned int iter;

    vector<double> dtList;
    vector<double> timeList;

    vector<ConservedVariable> convergenceHistory;

    long long computationTime;

// ====================================================================================================================
// ====================================================================================================================
//                      Methods
// ====================================================================================================================
// ====================================================================================================================
public:
	GKSSolver();

    GKSSolver(Parameters param, FluidParameter fluidParam);

	~GKSSolver();

    // ========================================================================
    //              Communication methods
    // ========================================================================

    void readMeshFromMeshObject( const GKSMesh& origin );

    void writeDataToMeshObject( const GKSMesh& target );

    // ========================================================================
    //              Simulation Control
    // ========================================================================

    void iterate();

    void timeStep();

    bool isConverged(ConservedVariable residual);

    void computeGlobalTimestep();

    void applyForcing(const idType id);

    void applyBoundaryCondition(const idType id);

    void computeFlux(const idType id);

    void updateCell(const idType id);

    // ========================================================================
    //              Flux computation subroutines
    // ========================================================================

    __declspec(noinline) PrimitiveVariable reconstructPrimPiecewiseConstant(const idType id);

    __declspec(noinline) ConservedVariable differentiateConsNormal(const idType id, double rho);

    void computeMicroSlope( const PrimitiveVariable& prim, const ConservedVariable& macroSlope, double* microSlope );

    void computeMoments(const PrimitiveVariable& prim, double * MomentU, double* MomentV, double * MomentXi, const int numberMoments);

    ConservedVariable computeTimeDerivative(double* MomentU, double* MomentV, double* MomentXi, double* a, double* b);

    ConservedVariable assembleFlux(double* MomentU, double* MomentV, double* MomentXi, double* a, double* b, double* A, double* timeCoefficients, PrimitiveVariable prim, double area, double tau);

    // ========================================================================
    //              Data Analysis
    // ========================================================================

    ConservedVariable getL2GlobalResidual();

    double getMaxVelocity();

    // ========================================================================
    //              Util
    // ========================================================================

    PrimitiveVariable cons2prim(ConservedVariable& cons);

    ConservedVariable prim2cons(PrimitiveVariable& prim);

    void local2global(const idType id, PrimitiveVariable& prim);
    void local2global(const idType id, ConservedVariable& cons);

    void global2local(const idType id, PrimitiveVariable& prim);
    void global2local(const idType id, ConservedVariable& cons);

    bool isGhostCell(const idType& id);

    idType findNeigborCellInDomain(const idType& id);

    ConservedVariable getData(idType id);
    PrimitiveVariable getPrim(idType id);

    void setData(idType id, ConservedVariable cons);
    void setData(idType id, PrimitiveVariable prim);
};

#endif