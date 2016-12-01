
#include "Types.h"
#include "BoundaryCondition.h"
#include <vector>
#include <array>
#include <fstream>

#ifndef MSHREADER_H
#define MSHREADER_H

using namespace std;

class mshReader
{
public:
    vector<Vec2> Nodes;
    
    vector< CellType >      Cell2Type;
    vector< array<idType,4> >  Cell2Node;
    vector< array<idType,4> >  Cell2Face;
    vector< idType >           Cell2BC;
                            
    vector< array<idType,2> >  Face2Cell;
    vector< array<idType,2> >  Face2Node;
    vector< idType >           Face2PhysicalName;
    vector< array<bool,2> > Face2CellAdd;

    vector<Vec2>   CellCenter;
    vector<double> CellVolume;
    vector<double> CellMinDx;
    vector< array<double,3> > CellLSCoeff;

    vector<Vec2>   FaceCenter;
    vector<Vec2>   FaceNormal;
    vector<double> FaceDistance;
    vector<double> FaceArea;

    vector< BoundaryCondition* > BCs;
    vector< string > BCNames;

    vector< idType > PhysicalNameIDs;
    vector< string > PhysicalNames;
    vector< idType > PhysicalNames2BCs;

public:
    mshReader();
    ~mshReader();

    bool readProblem( string problemName );

    bool readBoundaryConditions( string filename );

    bool readMsh(string filename);

    bool readPhysicalNames(ifstream& file);

    bool readNodes(ifstream& file);

    bool readElements(ifstream& file);

    bool newFace(string buffer);

    bool newCell(string buffer, CellType type);

    bool linkExistingFaces(array<idType,4> tmpCell2Node, array<idType,4>& tmpCell2Face, idType CellID, CellType type);

    bool createMissingFaces(array<idType,4> tmpCell2Node, array<idType,4>& tmpCell2Face, idType CellID, CellType type);

    void computeCellGeometry();

    void computeCellLSCoeff();

    void computeFaceGeometry();

    void computeCellMinDx();

    idType findPeriodicInterface( idType face );

    void createGhostCells();

    bool matchPhysicalNamesWithBCs();

    bool NodeCheck();

    bool FaceCheck();

    // ================================ Util ==================================

    template <typename T>
    idType findIndex(vector<T> _vector, T _value);

    double normalDistanceFace2Cell( idType face, idType cell );

    double normalDistanceFace2Point( idType face, Vec2 point );

    idType getNeighborCell( idType face, idType askingCell );
};

#endif

template<typename T>
inline idType mshReader::findIndex(vector<T> _vector, T _value)
{
    return find( _vector.begin(), _vector.end(), _value ) - _vector.begin();
}
