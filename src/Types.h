
#ifndef TYPES_H
#define TYPES_H

struct PrimitiveVariable
{
	double rho;
	double U;
	double V;
	double L;   // Lambda = m/(2kT) = 1/(2RT)
};

struct ConservedVariable
{
	double rho;
	double rhoU;
	double rhoV;
    double rhoE;
};

struct float2
{
	double x;
	double y;

    float2(){x = 0.0; y = 0.0;}

    float2(double xArg, double yArg)
    {
        x = xArg;
        y = yArg;
    }
};

struct FluidParameter
{
    int K;   
    double nu;  // viscosity
    double R;   // spez gasconstant
    float2 Force;
    float2 BoussinesqForce;
    double rhoReference;
};

struct Parameters
{
	unsigned int numberOfIterations;
    unsigned int outputInterval;
    unsigned int outputIntervalVTK;

    double convergenceCriterium[4];
    
    bool resOutput;
    bool fluxOutput;

    double L;
	double CFL;

    bool verbose;
};

enum InterfaceType
{
    compressible,
    incompressible
};

enum BoundaryConditionType
{
    wall,
    isothermalWall,
    periodic,
    periodicGhost
};

#endif