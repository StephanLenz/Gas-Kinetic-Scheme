// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      Cell.h
//
//      Function:
//          Storage of cell data
//          Implementation of cell related computations
//              Update of Conserved Variables
//              Forcing
//              Timestep computations
//
// ============================================================================

#ifndef RECTCELL2D_H
#define RECTCELL2D_H

//#include "Interface.h"
class Interface;
#include "Types.h"
#include "BoundaryCondition.h"
#include <string>
#include <array>

using namespace std;

class Cell
{
// ====================================================================================================================
// ====================================================================================================================
//                      Attributes
// ====================================================================================================================
// ====================================================================================================================
private:

    static unsigned long int numberOfCells;

    unsigned long int ID;

    InterfaceType interfaceType;        // compressible or incompressible

    FluidParameter fluidParam;

    // ========================================================================
    //          Connectivity as Pointers
    // ========================================================================
    Node* nodes[4];
    int nInterfaces;
	Interface* InterfaceList[4];
    BoundaryCondition* BoundaryContitionPointer;

    // ========================================================================
    //          Geometrical Information
    // ========================================================================
	Node center;
    double volume;
    double minDx;

    // ========================================================================
    //          Data
    // ========================================================================
	double cons[4];
    double cons_old[4];

    ConservedVariable updateVal;
    ConservedVariable residual;

    ConservedVariable gradientX;
    ConservedVariable gradientY;

    double r11, r12, r22;
    
// ====================================================================================================================
// ====================================================================================================================
//                      Methods
// ====================================================================================================================
// ====================================================================================================================
public:

	Cell();

	Cell(InterfaceType interfacetype, Node** nodes, BoundaryCondition* BC, FluidParameter fluidParam);

	~Cell();

    // ========================================================================
    //              Computation Methods
    // ========================================================================

    double getLocalTimestep();

    void applyForcing(double dt);

    void applyBoundaryCondition();
    
    void computeCons(PrimitiveVariable prim);

    void computeLeastSquareGradients();

	void update(double dt);

    // ========================================================================
    //              Initialization Methods
    // ========================================================================
    
    void addInterface(Interface* newInterface);

    void computeMinDx();

    void computeLeastSquareCoefficients();

	void setValues(double rho, double u, double v, double T);

    void setCons( ConservedVariable cons);
    
    // ========================================================================
    //              get Methods
    // ========================================================================

    // ============ get Connectivity ==========================================

    unsigned long int getID();

    Node getNode(int i);

    Interface* getInterface(int i);

    Cell* getNeighborCell(int i);

    Cell* findNeighborInDomain();

    Node getConnectivity(int i);

    bool isGhostCell();

    BoundaryCondition* getBoundaryConditionPointer();
    
    // ============ get Geometry ==============================================

    double getVolume();

	Node getCenter();

    double getMinDx();

    array<double,3> getLSCoeff();

    double distance(Node point);
    
    // ============ get Data ==================================================

    PrimitiveVariable getPrim();

    ConservedVariable getCons();

    ConservedVariable getConsOld();

    ConservedVariable getLocalResidual();

    ConservedVariable getUpdate();

    ConservedVariable getGradientX();

    ConservedVariable getGradientY();

    // ========================================================================
    //              output Methods
    // ========================================================================

	string toString();

    string valuesToString();

	string writeNodes();
};

#endif