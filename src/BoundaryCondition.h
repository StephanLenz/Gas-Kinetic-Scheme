
#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

class BoundaryCondition
{
private:
    short int type[4];
    double value[4];
public:
    BoundaryCondition();
    BoundaryCondition(  int rhoType, int UType, int VType, int TType,
                        double rho, double U, double V, double T);
    ~BoundaryCondition();

    short int getType(short int i);

    double getValue(short int i);
};

#endif

