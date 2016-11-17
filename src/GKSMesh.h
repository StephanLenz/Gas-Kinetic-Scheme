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

#ifndef GKSMESH_H
#define GKSMESH_H

#include "Cell.h"
#include "Interface.h"
#include "BoundaryCondition.h"
#include "InterfaceBC.h"
#include "Types.h"
#include <vector>

using namespace std;

class GKSMesh
{
// ====================================================================================================================
// ====================================================================================================================
//                      Attributes
// ====================================================================================================================
// ====================================================================================================================
public:

    // ========================================================================
    //              Mesh definition
    // ========================================================================
    vector<Node*>     NodeList;           // List of Pointers to nodes
	vector<Cell*>		CellList;           // List of Pointers to cells
	vector<Interface*>	InterfaceList;      // List of Pointers to interfaces

    vector<BoundaryCondition*> BoundaryConditionList;
    vector<InterfaceBC*> InterfaceBoundaryConditionsList;

	double lengthX;		// dimension of the domain in x direction
	double lengthY;     // dimension of the domain in y direction
    
    // ========================================================================
    //              Fluid Parameters
    // ========================================================================
    FluidParameter fluidParam;
    
    // ========================================================================
    //              Simulation Parameters/Data
    // ========================================================================
    Parameters param;
    double dt;
    vector<double> dtList;
    double time;
    vector<double> timeList;

    unsigned int iter;

    vector<ConservedVariable> convergenceHistory;

    long long computationTime;

// ====================================================================================================================
// ====================================================================================================================
//                      Methods
// ====================================================================================================================
// ====================================================================================================================
public:
	GKSMesh();

    GKSMesh(Parameters param, FluidParameter fluidParam);

	~GKSMesh();

    void generateRectMeshGraded(InterfaceType type, double lengthX, double lengthY, int nx, int ny, double gradingX, double gradingY);

    void addBoundaryCondition(  BoundaryConditionType type, double rho, double U, double V, double T);

    void addInterfaceBoundaryCondition(double wallVelocity);

    // ========================================================================
    //              Initialization methods
    // ========================================================================

	void initMeshConstant(double rho, double u, double v, double T);

    void initMeshLinear(double * rho, double * u, double * v, double * lambda);

    void initMeshLinearHorizontal(double * rho, double * u, double * v, double * lambda);

    void initMeshParabularVelocity(double rho, double u, double v, double T);

    void initMeshSineVelocity(double rho, double u, double v, double T);

    void initMeshStepVelocity(double rho, double u, double v, double T);

    void initMeshAtmospheric(double rho, double u, double v, double lambda, double g);

    // ========================================================================
    //              Simulation Control
    // ========================================================================

    void iterate();

    void timeStep();

    void computeGlobalTimestep();

    void applyForcing();

    void applyBoundaryCondition();

    void computeLeastSquareGradients();

    void computeFluxes();

    void updateCells();

    // ========================================================================
    //              Data Analysis
    // ========================================================================

    ConservedVariable getMaxGlobalResidual();
    ConservedVariable getL2GlobalResidual();

    double getMaxVelocity();

    double getMaxRe();
    double getMaxMa();

    bool isConverged(ConservedVariable residual);

    // ========================================================================
    //              methods for file output
    // ========================================================================

    void writeOutputFiles();

    void writeOverviewFile(string filename);

    void writeTimeSteps(string filename);

    void writeTime(string filename);

    void writeConvergenceHistory(string filename);

	void writeVTKFile(string filename, bool data = true, bool BC = false);

    void writeVTKFileFlux(string filename, bool data = true, bool BC = false);

    void writeResultFields(string filename);

    void writeResultBoundaryFluxes(string filename);

    void writeGambitNeutralFile(string filename);

    // ========================================================================
    //              private methods for VTK output
    // ========================================================================
private:

    void writeCellGeometry(ofstream& file);

    void writeInterfaceGeometry(ofstream& file);

    void writeCellData(ofstream& file);

    void writeInterfaceData(ofstream& file);

};

#endif