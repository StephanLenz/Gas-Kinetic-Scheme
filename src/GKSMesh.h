
#ifndef GKSMESH_H
#define GKSMESH_H

#include "Cell.h"
#include "Interface.h"
#include "BoundaryCondition.h"
#include "InterfaceBC.h"
#include "Types.h"
#include <vector>
#include <string>
#include <chrono>

using namespace std;

using namespace std;

class GKSMesh
{
private:
	vector<Cell*>		CellList;
	vector<Interface*>	InterfaceList;

    vector<BoundaryCondition*> BoundaryConditionList;
    vector<InterfaceBC*> InterfaceBoundaryConditionsList;

    Parameters param;
    FluidParameter fluidParam;

	double lengthX;		
	double lengthY;

    double dt;
    vector<double> dtList;
    double time;
    unsigned int iter;

    vector<ConservedVariable> convergenceHistory;

    long long computationTime;

public:
	GKSMesh();

    GKSMesh(Parameters param, FluidParameter fluidParam);

	~GKSMesh();

	void generateRectMesh(InterfaceType type, double lengthX, double lengthY, int nx, int ny);

    void generateRectMeshPeriodic(InterfaceType type, double lengthX, double lengthY, int nx, int ny);

    void generateRectMeshPeriodicVertical(InterfaceType type, double lengthX, double lengthY, int nx, int ny);

    void generateRectMeshInterfaceBCs(InterfaceType type, double lengthX, double lengthY, int nx, int ny);

    void generateRectMeshPeriodicInterfaceBCs(InterfaceType type, double lengthX, double lengthY, int nx, int ny);

	void initMeshConstant(double rho, double u, double v, double T);

	void initMeshLinearTemperature(double rho, double u, double v, double * T);

    void initMeshLinear(double * rho, double * u, double * v, double * lambda);

    void initMeshLinearDensity(double* rho, double u, double v, double T);

    void addBoundaryCondition(  int rhoType, int UType, int VType, int TType,
                                double rho, double U, double V, double T);

    void addInterfaceBoundaryCondition(double wallVelocity);

    void applyBoundaryCondition();

    void computeGlobalTimestep();

    ConservedVariable getMaxGlobalResidual();
    ConservedVariable getL2GlobalResidual();

    void timeStep();

    void iterate();

    double getMaxVelocity();

    double getMaxRe();
    double getMaxMa();

    bool isConverged(ConservedVariable residual);

	string toString();
    string cellValuesToString();

    void writeOverviewFile(string filename);

	void writeVTKFile(string filename, bool data = true, bool BC = false);

    void writeVTKFileFlux(string filename, bool data = true, bool BC = false);

    void writeTimeSteps(string filename);

    void writeVelocityProfile(string filename, double x);
    void writePressureGradientProfile(string filename, double x);

    void writeMeshAsText(string filename);

    void writeVelocityU(string filename);

    void writeVelocityV(string filename);

    void writeConvergenceHistory(string filename);

private:

    void writeCellGeometry(ofstream& file);

    void writeInterfaceGeometry(ofstream& file);

    void writeCellData(ofstream& file);

    void writeInterfaceData(ofstream& file);
};

#endif