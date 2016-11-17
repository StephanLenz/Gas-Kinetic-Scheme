
#include "Types.h"
#include <vector>
#include <array>
#include <fstream>

#ifndef MSHREADER_H
#define MSHREADER_H

using namespace std;

class mshReader
{
    //struct Vec2 { 
    //    double x; 
    //    double y; 
    //    Vec2(double x, double y):x(x),y(y){} 
    //    Vec2():x(0.0),y(0.0){}
    //};

    //enum BCtype { periodicLeft, periodicRight, wall};

public:

    vector<Vec2> Nodes;
    
    vector< array<int,3> >  Cell2Node;
    vector< array<int,3> >  Cell2Face;
    vector< int >           Cell2BC;
                            
    vector< array<int,2> >  Face2Cell;
    vector< array<int,2> >  Face2Node;
    vector< int >           Face2BC;
    vector< array<bool,2> > Face2CellAdd;

    vector<Vec2>   CellCenter;
    vector<double> CellVolume;
    vector<double> CellMinDx;

    vector<Vec2>   FaceCenter;
    vector<Vec2>   FaceNormal;
    vector<double> FaceDistance;
    vector<double> FaceArea;

    vector< int >    BCIDs;
    vector< BoundaryConditionType > BCs;

public:
    mshReader();
    ~mshReader();

    bool readMsh(string filename);

    bool readPhysicalNames(ifstream& file);

    bool readNodes(ifstream& file);

    bool readElements(ifstream& file);

    bool newFace(string buffer);

    bool newCell(string buffer);

    bool linkExistingFaces(array<int,3> tmpCell2Node, array<int,3>& tmpCell2Face, int CellID);

    bool createMissingFaces(array<int,3> tmpCell2Node, array<int,3>& tmpCell2Face, int CellID);

    void computeCellGeometry();

    void computeFaceGeometry();

    void computeCellMinDx();

    bool findPeriodicCells();

    void createGhostCells();

    template <typename T>
    int findIndex(vector<T> _vector, T _value);

    double normalDistanceFace2Cell( int face, int cell );

    double normalDistanceFace2Point( int face, Vec2 point );
};

#endif

template<typename T>
inline int mshReader::findIndex(vector<T> _vector, T _value)
{
    return find( _vector.begin(), _vector.end(), _value ) - _vector.begin();
}
