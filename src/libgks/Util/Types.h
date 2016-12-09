#include <string>

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

    ConservedVariable() : rho(0.0), rhoU(0.0), rhoV(0.0), rhoE(0.0) {}
};

struct Node
{
	double x;
	double y;
    long int ID;

    Node(){x = 0.0; y = 0.0;}

    Node(double xArg, double yArg)
    {
        x = xArg;
        y = yArg;
    }
};


struct Vec2
{
	double x;
	double y;

    Vec2(){x = 0.0; y = 0.0;}

    Vec2(double xArg, double yArg)
    {
        x = xArg;
        y = yArg;
    }

    Vec2& operator=( const Node& that )
    {
        this->x = that.x;
        this->y = that.y;
        return *this;
    }

};

struct FluidParameter
{
    double K;   
    double nu;  // viscosity
    double R;   // spez gasconstant
    Node Force;
    Node relativeForce;
    PrimitiveVariable referencePrim;
    double Pr;
};

struct Parameters
{
    std::string simulationName;

	unsigned long int numberOfIterations;
    unsigned long int outputInterval;
    unsigned long int outputIntervalVTK;

    double maxTime;

    ConservedVariable convergenceCriterium;

	double CFL;
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
    periodicGhost,
    inlet,
    outlet,
    none
};

enum CellType
{
    tri,
    quad
};

typedef long int idType;

#endif
