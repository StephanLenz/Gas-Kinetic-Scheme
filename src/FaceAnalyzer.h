#include "Types.h"
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
};

class DragAndLiftCalculator : public FaceAnalyzer
{
private:
    double drag;
    double lift;
public:
    DragAndLiftCalculator( );
    virtual void analyze( GKSSolver& solver );
    virtual void print();
};

#endif

