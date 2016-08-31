
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
    vector<float2*>     NodeList;
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

    void generateRectMeshGraded(InterfaceType type, double lengthX, double lengthY, int nx, int ny, double gradingX, double gradingY);

    void generateMiniPatchMesh();

	void initMeshConstant(double rho, double u, double v, double T);

	void initMeshLinearTemperature(double rho, double u, double v, double * T);

    void initMeshLinear(double * rho, double * u, double * v, double * lambda);

    void initMeshLinearHorizontal(double * rho, double * u, double * v, double * lambda);

    void initMeshLinearDensity(double* rho, double u, double v, double T);

    void initMeshParabularVelocity(double rho, double u, double v, double T);

    void initMeshSineVelocity(double rho, double u, double v, double T);

    void initMeshAtmospheric(double rho, double u, double v, double lambda, double g);


    void addBoundaryCondition(  BoundaryConditionType type,
                                double rho, double U, double V, double T);

    void addInterfaceBoundaryCondition(double wallVelocity);

    void applyBoundaryCondition();

    void applyForcing();

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
    void writeResultFields(string filename);
    void writeTemperatureProfile(string filename, double x);

    void writeMeshAsText(string filename);

    void writeVelocityU(string filename);

    void writeVelocityV(string filename);

    void writeTemperature(string filename);

    void writeDensity(string filename);

    void writeConvergenceHistory(string filename);

private:

    void writeCellGeometry(ofstream& file);

    void writeInterfaceGeometry(ofstream& file);

    void writeCellData(ofstream& file);

    void writeInterfaceData(ofstream& file);
};

#endif