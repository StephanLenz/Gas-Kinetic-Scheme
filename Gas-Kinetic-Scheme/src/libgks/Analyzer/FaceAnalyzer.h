#include "libgks/Util/Types.h"
class GKSSolver;
#include <string>
#include <vector>
#ifndef FACEANALYZER_H
#define FACEANALYZER_H

using namespace std;

class FaceAnalyzer
{
protected:
    vector<idType> Faces;

public:
    void addFace(idType id);
    virtual void analyze( GKSSolver& solver ) = 0;
    virtual void print( ) = 0;
    virtual void write( string filename, double t ) = 0;
};

class DragAndLiftCalculator : public FaceAnalyzer
{
private:
    double drag;
    double lift;
    double U;
    double D;
public:
    DragAndLiftCalculator( double U, double D );
    virtual void analyze( GKSSolver& solver );
    virtual void print();
    virtual void write( string filename, double t );
private:
    ConservedVariable computeTimeAveragedFlux( GKSSolver& solver, idType face );
};

#endif

