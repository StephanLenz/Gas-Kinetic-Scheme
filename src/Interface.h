
#ifndef INTERFACE_H
#define INTERFACE_H

#include "Cell.h"
#include "InterfaceBC.h"
//class Cell;
#include <string>

using namespace std;

class Interface
{
protected:
    static unsigned long int numberOfCells;
    unsigned long int ID;

public:
	Cell* negCell;
	Cell* posCell;
protected:
    float2* nodes[2];
    float2 center;
    float2 normal;

    double area;

    double posDistance;
    double negDistance;

    BoundaryCondition* BoundaryConditionPointer;

    FluidParameter fluidParam;

    double timeIntegratedFlux[4];
    double timeIntegratedFlux_1[4];
    double timeIntegratedFlux_2[4];
    double timeIntegratedFlux_3[4];
    double FluxDensity[4];

    static int interpolationOrder;

public:
	Interface();
	Interface(Cell* negCell, Cell* posCell, bool negAdd, bool posAdd,
              float2** nodes, FluidParameter fluidParam, BoundaryCondition* BC, double periodicLengthX, double periodicLengthY);
	~Interface();

    static Interface* createInterface(InterfaceType type, Cell* negCell, Cell* posCell, bool negAdd, bool posAdd,
                                      float2** nodes, FluidParameter fluidParam, BoundaryCondition* BC, double periodicLengthX, double periodicLengthY);


	virtual void computeFlux(double dt);

    virtual void computeInternalFlux(double dt);

    Cell* getNeigborCell(Cell* askingCell);
    Cell* getCellInDomain();

    ConservedVariable getTimeIntegratedFlux();
    ConservedVariable getTimeIntegratedFlux_1();
    ConservedVariable getTimeIntegratedFlux_2();
    ConservedVariable getTimeIntegratedFlux_3();
    ConservedVariable getFluxDensity();
    double getFluxSign(Cell* askingCell);

    bool isGhostInterface();
    bool isBoundaryInterface();

    void addCell(Cell* that);
    Cell* getPeriodicCell();

    float2* getNode(int i);
    float2 getNormal();
    float2 getCenter();
    float2 getScaledNormal();
    double getArea();

    BoundaryCondition* getBoundaryCondition();
    float2 getPosConnectivity();
    float2 getNegConnectivity();


    double distance(float2 point);

	string toString();

    string writeCenter();

    static void setInterpolationOrder(int arg );

    static int getInterpolationOrder();

protected:

    void interpolatePrim(double* prim);

    void reconstructPrimPiecewiseConstant(double* prim);

    void reconstructPrimPiecewiseLinear(double* prim);

    void differentiateConsNormal(double* normalGradCons, double* prim);

    void differentiateConsNormalThreePoint(double* normalGradCons, double* prim);

    void differentiateConsLeastSquare(double* normalGradCons, double* tangentialGradCons, double* prim);

    virtual void computeTimeDerivative(double* prim, double* MomentU, double* MomentV, double* MomentXi,
                                       double* a, double* b, double * timeGrad) = 0;

    virtual void assembleFlux(double* MomentU, double* MomentV, double* MomentXi, 
                              double* a, double* b, double* A, double* timeCoefficients,
                              double* prim, double tau) = 0;

    void transformGlobal2Local(double* vec);
    void transformLocal2Global(double * vec);

    PrimitiveVariable cons2Prim(ConservedVariable cons);

    virtual void computeMicroSlope(double* prim, double* macroSlope, double* microSlope) = 0;
    void computeMoments(double* prim, double* MomentU, double* MomentV, double* MomentXi, int numberMoments);
};

#endif