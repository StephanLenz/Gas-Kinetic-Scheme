#include "libgks/Util/Types.h"
class mshReader;
class GKSSolver;
#include <string>
#include <vector>
#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

using namespace std;

class BoundaryCondition
{
protected:
    vector<idType> Cell;
    vector<idType> NeighborCell;

public:
    virtual ~BoundaryCondition();

    void addCell(idType id);
    virtual void findNeighborCells( mshReader& reader );

    void setGhostCells( GKSSolver& solver );

    virtual void setGhostCell(GKSSolver& solver, idType cell) = 0;
    virtual void setGradientGhostCells(GKSSolver& solver) = 0;
};

class bcWall : public BoundaryCondition
{
private:
    double U;
    double V;
public:
    virtual ~bcWall();
    bcWall( double U, double V );
    virtual void setGhostCell(GKSSolver& solver, idType cell);
    virtual void setGradientGhostCells(GKSSolver& solver);
};

class bcIsothermalWall : public BoundaryCondition
{
private:
    double U;
    double V;
    double L;
public:
    virtual ~bcIsothermalWall();
    bcIsothermalWall( double U, double V, double L );
    virtual void setGhostCell(GKSSolver& solver, idType cell);
    virtual void setGradientGhostCells(GKSSolver& solver);
};

class bcPeriodicGhost : public BoundaryCondition
{
public:
    virtual ~bcPeriodicGhost();
    bcPeriodicGhost( );
    virtual void findNeighborCells( mshReader& reader );
    virtual void setGhostCell(GKSSolver& solver, idType cell);
    virtual void setGradientGhostCells(GKSSolver& solver);
};

class bcInflowParabolic : public BoundaryCondition
{
private:
    double U;
    double V;
    double L;
    Vec2 p0;
    Vec2 p1;
public:
    virtual ~bcInflowParabolic();
    bcInflowParabolic( double U, double V, double L, Vec2 p0, Vec2 p1 );
    virtual void setGhostCell(GKSSolver& solver, idType cell);
    virtual void setGradientGhostCells(GKSSolver& solver);
};

class bcInflowUniform : public BoundaryCondition
{
private:
    double U;
    double V;
    double L;
public:
    virtual ~bcInflowUniform();
    bcInflowUniform( double U, double V, double L );
    virtual void setGhostCell(GKSSolver& solver, idType cell);
    virtual void setGradientGhostCells(GKSSolver& solver);
};

class bcOutflow : public BoundaryCondition
{
private:
    double p;
public:
    virtual ~bcOutflow();
    bcOutflow( double p );
    virtual void setGhostCell(GKSSolver& solver, idType cell);
    virtual void setGradientGhostCells(GKSSolver& solver);
};

#endif
