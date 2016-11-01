// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      Interface.h
//
//          This class in an abstract baseclass for compressible and
//          incompressible Interfaces. It implements all methods that do not
//          depend on the type of interface.
//
//      Function:
//          Storage of Interface data
//          Implementation of Interface related computations
//              Flux computation 
//
// ============================================================================

#ifndef INTERFACE_H
#define INTERFACE_H

#include "Cell.h"
#include "InterfaceBC.h"
//class Cell;
#include <string>

using namespace std;

class Interface
{
// ====================================================================================================================
// ====================================================================================================================
//                      Attributes
// ====================================================================================================================
// ====================================================================================================================
protected:
    static unsigned long int numberOfCells;
    unsigned long int ID;

    static int interpolationOrder;

    FluidParameter fluidParam;
    
    // ========================================================================
    //          Connectivity as Pointers
    // ========================================================================
    Node* nodes[2];

	Cell* negCell;
	Cell* posCell;

    BoundaryCondition* BoundaryConditionPointer;

    // ========================================================================
    //          Geometrical Information
    // ========================================================================
    Node center;
    Node normal;

    double area;

    double posDistance;
    double negDistance;

    // ========================================================================
    //          Data
    // ========================================================================
    double timeIntegratedFlux[4];
    double timeIntegratedFlux_1[4];
    double timeIntegratedFlux_2[4];
    double timeIntegratedFlux_3[4];
    double FluxDensity[4];
    
// ====================================================================================================================
// ====================================================================================================================
//                      Methods
// ====================================================================================================================
// ====================================================================================================================
public:
	Interface();

	Interface(Cell* negCell, Cell* posCell, bool negAdd, bool posAdd,
              Node** nodes, FluidParameter fluidParam, BoundaryCondition* BC,
              double periodicLengthX, double periodicLengthY);

	~Interface();

    // ========================================================================
    //              Interface Factory
    // ========================================================================
    static Interface* createInterface(InterfaceType type, Cell* negCell, Cell* posCell, bool negAdd, bool posAdd,
                                      Node** nodes, FluidParameter fluidParam, BoundaryCondition* BC,
                                      double periodicLengthX, double periodicLengthY);
    
    // ========================================================================
    //              Initialization Methods
    // ========================================================================

    void addCell(Cell* that);

    static void setInterpolationOrder(int arg );

    static int getInterpolationOrder();

    // ========================================================================
    //              Computation Methods
    // ========================================================================
public:
	virtual void computeFlux(double dt);

protected:
    // ============ Reconstruction of primitive Variables =====================

    void interpolatePrim(double* prim);

    __declspec(noinline) void reconstructPrimPiecewiseConstant(double* prim);

    void reconstructPrimPiecewiseLinear(double* prim);

    // ============ Differentiation of conserved Variables ====================

    __declspec(noinline) void differentiateConsNormal(double* normalGradCons, double* prim);

    void differentiateConsNormalThreePoint(double* normalGradCons, double* prim);

    void differentiateConsLeastSquare(double* normalGradCons, double* tangentialGradCons, double* prim);

    // ============ Flux Computation ==========================================

    virtual void computeInternalFlux(double dt);

    void computeMoments(double* prim, double* MomentU, double* MomentV, double* MomentXi, int numberMoments);

    virtual void computeMicroSlope(double* prim, double* macroSlope, double* microSlope) = 0;

    virtual void computeTimeDerivative(double* prim, double* MomentU, double* MomentV, double* MomentXi,
                                       double* a, double* b, double * timeGrad) = 0;

    virtual void assembleFlux(double* MomentU, double* MomentV, double* MomentXi, 
                              double* a, double* b, double* A, double* timeCoefficients,
                              double* prim, double tau) = 0;

    // ============ Data Manipulation =========================================

    void transformGlobal2Local(double* vec);

    void transformLocal2Global(double * vec);

    PrimitiveVariable cons2Prim(ConservedVariable cons);

    // ========================================================================
    //              get Methods
    // ========================================================================
public:
    // ============ get Connectivity ==========================================

    idType getID();

    Node* getNode(int i);

    Cell* getCell( idType i );

    Cell* getNeigborCell(Cell* askingCell);

    Cell* getCellInDomain();

    Cell* getPeriodicCell();

    Node getPosConnectivity();

    Node getNegConnectivity();

    BoundaryCondition* getBoundaryCondition();

    bool isGhostInterface();

    bool isBoundaryInterface();

    // ============ get Geometry ==============================================

    Node getNormal();

    Node getCenter();

    Node getScaledNormal();

    double getArea();

    double getDistance2CellCenter(idType i);

    double distance(Node point);
    
    // ============ get Data ==================================================

    ConservedVariable getTimeIntegratedFlux();

    ConservedVariable getTimeIntegratedFlux_1();

    ConservedVariable getTimeIntegratedFlux_2();

    ConservedVariable getTimeIntegratedFlux_3();

    ConservedVariable getFluxDensity();

    double getFluxSign(Cell* askingCell);

    // ========================================================================
    //              output Methods
    // ========================================================================

	string toString();

    string writeCenter();

protected:
};

#endif