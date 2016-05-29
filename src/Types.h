
#ifndef TYPES_H
#define TYPES_H

struct PrimaryVariable
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
};

struct FluidParameter
{
    int K;   
    double nu;  // viscosity
    double R;   // spez gasconstant
    float2 Force;
};

struct Parameters
{
	unsigned int numberOfIterations;
    unsigned int outputInterval;

	double CFL;

    bool verbose;
};

#endif