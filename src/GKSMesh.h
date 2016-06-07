
#ifndef GKSMESH_H
#define GKSMESH_H

#include "Cell.h"
#include "Interface.h"
#include "BoundaryCondition.h"
#include "Types.h"
#include <vector>
#include <string>

using namespace std;

using namespace std;

class GKSMesh
{
private:
	vector<Cell*>		CellList;
	vector<Interface*>	InterfaceList;
    vector<BoundaryCondition*> BoundaryConditionList;

    Parameters param;
    FluidParameter fluidParam;

	double lengthX;		
	double lengthY;

    double dt;
    vector<double> dtList;
    unsigned int iter;

public:
	GKSMesh();

    GKSMesh(Parameters param, FluidParameter fluidParam);

	~GKSMesh();

	void generateRectMesh(double lengthX, double lengthY, int nx, int ny);

    void generateRectMeshPeriodic(double lengthX, double lengthY, int nx, int ny);

	void initMeshConstant(double rho, double u, double v, double T);

	void initMeshLinearTemperature(double rho, double u, double v, double * T);

    void initMeshLinearDensity(double* rho, double u, double v, double T);

    void addBoundaryCondition(  int rhoType, int UType, int VType, int TType,
                                double rho, double U, double V, double T);

    void applyBoundaryCondition();

    void computeGlobalTimestep();

    void timeStep();

    void iterate();

	string toString();
    string cellValuesToString();

	void writeVTKFile(string filename, bool data = true, bool BC = false);

    void writeVTKFileFlux(string filename, bool data = true, bool BC = false);

    void writeTimeSteps(string filename);

    void writeVelocityProfile(string filename, double x);

    void writeMeshAsText(string filename);

private:

    void writeCellGeometry(ofstream& file);

    void writeInterfaceGeometry(ofstream& file);

    void writeCellData(ofstream& file);

    void writeInterfaceData(ofstream& file);
};

#endif