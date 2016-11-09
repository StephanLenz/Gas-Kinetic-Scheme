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

#ifndef GKSSolver_H
#define GKSSolver_H

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

class GKSSolver
{
// ====================================================================================================================
// ====================================================================================================================
//                      Attributes
// ====================================================================================================================
// ====================================================================================================================
protected:
    idType numberOfNodes;
    idType numberOfCells;
    idType numberOfInterfaces;

    vector<BoundaryCondition> BoundaryConditionList;

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

    virtual void readMeshFromMeshObject( const GKSMesh& origin ) = 0;

    virtual void writeDataToMeshObject( const GKSMesh& target ) = 0;

    // ========================================================================
    //              Simulation Control
    // ========================================================================

    void iterate();

    void timeStep();

    bool isConverged(ConservedVariable residual);

    void computeGlobalTimestep();

    void applyForcing(const idType id);

    void computeFlux(const idType id);

    virtual void updateCell(const idType id) = 0;

    // ========================================================================
    //              Flux computation subroutines
    // ========================================================================

    virtual void storeDataOld(idType id) = 0;

    PrimitiveVariable reconstructPrimPiecewiseConstant(const idType id);

    ConservedVariable differentiateConsNormal(const idType id, double rho);

    void computeMicroSlope( const PrimitiveVariable& prim, const ConservedVariable& macroSlope, double* microSlope );

    void computeMoments(const PrimitiveVariable& prim, double * MomentU, double* MomentV, double * MomentXi, const int numberMoments);

    ConservedVariable computeTimeDerivative(double* MomentU, double* MomentV, double* MomentXi, double* a, double* b);

    ConservedVariable assembleFlux(double* MomentU, double* MomentV, double* MomentXi, double* a, double* b, double* A, double* timeCoefficients, PrimitiveVariable prim, double area, double tau);

    virtual void applyFlux(idType id, ConservedVariable flux) = 0;

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

    virtual bool isGhostCell(const idType& id) = 0;

    virtual ConservedVariable getCellData(idType id) = 0;
    virtual ConservedVariable getCellDataOld(idType id) = 0;

    virtual double getCellMinDx(idType id) = 0;

    virtual double getInterfaceArea(idType id) = 0;
    virtual double getInterfaceDistance(idType id) = 0;

    virtual idType getPosCell(idType id) = 0;
    virtual idType getNegCell(idType id) = 0;

    virtual Vec2 getInterfaceNormal(idType id) = 0;

    virtual void setData(idType id, ConservedVariable cons) = 0;
    virtual void setData(idType id, PrimitiveVariable prim);

    virtual PrimitiveVariable  getPrim(idType id);
};

#endif