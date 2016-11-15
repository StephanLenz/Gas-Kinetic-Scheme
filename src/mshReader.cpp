#include "mshReader.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

mshReader::mshReader()
{
}


mshReader::~mshReader()
{
}

bool mshReader::readMsh(string filename)
{
    ifstream file;
    file.open(filename);

    string buffer;

    getline(file, buffer);
    getline(file, buffer);
    getline(file, buffer);

    if( ! this->readPhysicalNames(file) ) return false;
    if( ! this->readNodes(file) )         return false;
    if( ! this->readElements(file) )      return false;
    
    if( ! this->findPeriodicCells() )     return false;

    return true;
}

bool mshReader::readPhysicalNames(ifstream& file)
{
    cout << "Read Boundary Conditions:" << endl;

    string buffer;
    getline(file, buffer);
    if( buffer.compare("$PhysicalNames") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }

    int nPhysicalNames;
    file >> nPhysicalNames;
    getline(file, buffer); // go to next line

    for( int i = 0; i < nPhysicalNames; ++i )
    {
        getline(file, buffer);
        stringstream bufferStream(buffer);
        
        int PhysicalNameType, ID;
        string type;

        bufferStream >> PhysicalNameType;
        if(PhysicalNameType != 1) continue;

        bufferStream >> ID >> type;

        this->BCIDs.push_back(ID);

        type = type.substr(1, type.length()-2 );

        if     ( type.compare("periodicLeft")  == 0 ) this->BCs.push_back(periodicLeft);
        else if( type.compare("periodicRight") == 0 ) this->BCs.push_back(periodicRight);
        else if( type.compare("wall")          == 0 ) this->BCs.push_back(wall);
        else { cout << "Error: Wrong BC" << endl; return false; }
    }

    getline(file, buffer);
    if( buffer.compare("$EndPhysicalNames") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }
    
    return true;
}

bool mshReader::readNodes(ifstream & file)
{
    cout << "Read Node Coordinates:" << endl;

    string buffer;
    getline(file, buffer);
    if( buffer.compare("$Nodes") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }

    int nNodes;
    file >> nNodes;
    getline(file, buffer); // go to next line

    for( int i = 0; i < nNodes; ++i )
    {
        getline(file, buffer);
        stringstream bufferStream(buffer);
        
        int ID;
        Vec2 tmpNode(0,0);

        bufferStream >> ID;
        if(ID != i+1) { cout << "Error: Invalid Node ID" << endl; return false; }

        bufferStream >> tmpNode.x >> tmpNode.y;

        this->Nodes.push_back(tmpNode);
    }

    getline(file, buffer);
    if( buffer.compare("$EndNodes") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }

    return true;
}

bool mshReader::readElements(ifstream & file)
{
    cout << "Read Cells and Boundary Interfaces:" << endl;

    string buffer;
    getline(file, buffer);
    if( buffer.compare("$Elements") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }

    int nElements;
    file >> nElements;
    getline(file, buffer); // go to next line

    for( int i = 0; i < nElements; ++i )
    {
        getline(file, buffer);
        stringstream bufferStream(buffer);
        
        int ID, type;

        bufferStream >> ID >> type;

        if (type == 1) if( ! newFace(buffer) ) return false;
        if (type == 2) if( ! newCell(buffer) ) return false;
    }

    getline(file, buffer);
    if( buffer.compare("$EndElements") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }

    return true;
}

bool mshReader::newFace(string buffer)
{
    stringstream bufferStream(buffer);

    int ID;
    int BCID;
    int N1, N2;
    int tmp;

    bufferStream >> ID >> tmp >> tmp >> BCID >> tmp >> N1 >> N2;
    
    for(int i = 0; i < this->Face2Node.size(); ++i)
        if (  ( this->Face2Node[i][0] == N1-1 && this->Face2Node[i][1] == N2-1 )
           || ( this->Face2Node[i][1] == N1-1 && this->Face2Node[i][0] == N2-1 ) )
        {
            this->Face2BC[i] = findIndex( this->BCIDs, BCID );
            return true;
        }

    array<int,2> tmpFace2Node;
    tmpFace2Node[0] = N1-1;
    tmpFace2Node[1] = N2-1;
    this->Face2Node.push_back(tmpFace2Node);

    array<int,2> tmpFace2Cell;
    tmpFace2Cell[0] = -1;
    tmpFace2Cell[1] = -1;
    this->Face2Cell.push_back(tmpFace2Cell);

    array<bool,2> tmpFace2CellAdd = {false, false};
    this->Face2CellAdd.push_back(tmpFace2CellAdd);
    
    this->Face2BC.push_back( findIndex( this->BCIDs, BCID ) );

    return true;
}

bool mshReader::newCell(string buffer)
{
    stringstream bufferStream(buffer);

    int ID;
    int BCID;
    int N1, N2, N3;
    int tmp;

    bufferStream >> ID >> tmp >> tmp >> tmp >> tmp >> N1 >> N2 >> N3;

    array<int,3> tmpCell2Node;
    tmpCell2Node[0] = N1-1;
    tmpCell2Node[1] = N2-1;
    tmpCell2Node[2] = N3-1;

    array<int, 3> tmpCell2Face;
    tmpCell2Face[0] = -1;
    tmpCell2Face[1] = -1;
    tmpCell2Face[2] = -1;

    this->linkExistingFaces( tmpCell2Node, tmpCell2Face, this->Cell2Node.size() );

    this->createMissingFaces( tmpCell2Node, tmpCell2Face, this->Cell2Node.size() );

    this->Cell2Node.push_back(tmpCell2Node);
    this->Cell2Face.push_back(tmpCell2Face);

    return true;
}

bool mshReader::linkExistingFaces(array<int, 3> tmpCell2Node, array<int, 3>& tmpCell2Face, int CellID)
{
    for( int i = 0; i < Face2Node.size(); ++i )
    {            
        if      ( Face2Node[i][0] == tmpCell2Node[0] && Face2Node[i][1] == tmpCell2Node[1] )
        {
            tmpCell2Face[0] = i;
            Face2Cell   [i][1] = CellID;
            Face2CellAdd[i][1] = true;
        }
        else if ( Face2Node[i][1] == tmpCell2Node[0] && Face2Node[i][0] == tmpCell2Node[1] )
        {
            tmpCell2Face[0] = i;
            Face2Cell   [i][0] = CellID;
            Face2CellAdd[i][0] = true;
        }   
        // ======
        if      ( Face2Node[i][0] == tmpCell2Node[1] && Face2Node[i][1] == tmpCell2Node[2] )
        {
            tmpCell2Face[1] = i;
            Face2Cell   [i][1] = CellID;
            Face2CellAdd[i][1] = true;
        }
        else if ( Face2Node[i][1] == tmpCell2Node[1] && Face2Node[i][0] == tmpCell2Node[2] )
        {
            tmpCell2Face[1] = i;
            Face2Cell   [i][0] = CellID;
            Face2CellAdd[i][0] = true;
        } 
        // ======
        if      ( Face2Node[i][0] == tmpCell2Node[2] && Face2Node[i][1] == tmpCell2Node[0] )
        {
            tmpCell2Face[2] = i;
            Face2Cell   [i][1] = CellID;
            Face2CellAdd[i][1] = true;
        }
        else if ( Face2Node[i][1] == tmpCell2Node[2] && Face2Node[i][0] == tmpCell2Node[0] )
        {
            tmpCell2Face[2] = i;
            Face2Cell   [i][0] = CellID;
            Face2CellAdd[i][0] = true;
        }
    }

    return true;
}

bool mshReader::createMissingFaces(array<int, 3> tmpCell2Node, array<int, 3>& tmpCell2Face, int CellID)
{
    if(tmpCell2Face[0] == -1)
    {
        array<int,2> tmpFace2Node;
        tmpFace2Node[0] = tmpCell2Node[0];
        tmpFace2Node[1] = tmpCell2Node[1];
        Face2Node.push_back(tmpFace2Node);

        array<int,2> tmpFace2Cell;
        tmpFace2Cell[0] = -1;
        tmpFace2Cell[1] = CellID;
        Face2Cell.push_back(tmpFace2Cell);

        array<bool,2> tmpFace2CellAdd = {true, true};
        this->Face2CellAdd.push_back(tmpFace2CellAdd);

        Face2BC.push_back(-1);

        tmpCell2Face[0] = Face2Cell.size()-1;
    }

    if(tmpCell2Face[1] == -1)
    {
        array<int,2> tmpFace2Node;
        tmpFace2Node[0] = tmpCell2Node[1];
        tmpFace2Node[1] = tmpCell2Node[2];
        Face2Node.push_back(tmpFace2Node);

        array<int,2> tmpFace2Cell;
        tmpFace2Cell[0] = -1;
        tmpFace2Cell[1] = CellID;
        Face2Cell.push_back(tmpFace2Cell);

        array<bool,2> tmpFace2CellAdd = {true, true};
        this->Face2CellAdd.push_back(tmpFace2CellAdd);

        Face2BC.push_back(-1);

        tmpCell2Face[1] = Face2Cell.size()-1;
    }

    if(tmpCell2Face[2] == -1)
    {
        array<int,2> tmpFace2Node;
        tmpFace2Node[0] = tmpCell2Node[2];
        tmpFace2Node[1] = tmpCell2Node[0];
        Face2Node.push_back(tmpFace2Node);

        array<int,2> tmpFace2Cell;
        tmpFace2Cell[0] = -1;
        tmpFace2Cell[1] = CellID;
        Face2Cell.push_back(tmpFace2Cell);

        array<bool,2> tmpFace2CellAdd = {true, true};
        this->Face2CellAdd.push_back(tmpFace2CellAdd);

        Face2BC.push_back(-1);

        tmpCell2Face[2] = Face2Cell.size()-1;
    }
    
    return true;
}

void mshReader::computeCellGeometry()
{
    this->CellCenter.resize( this->Cell2Node.size() );
    this->CellVolume.resize( this->Cell2Node.size() );

    for( int cell = 0; cell < this->Cell2Node.size(); ++cell )
    {
        this->CellCenter[cell].x =  (this->Nodes[ Cell2Node[cell][0] ].x + this->Nodes[ Cell2Node[cell][1] ].x + this->Nodes[ Cell2Node[cell][2] ].x) / 3.0;
        this->CellCenter[cell].y =  (this->Nodes[ Cell2Node[cell][0] ].y + this->Nodes[ Cell2Node[cell][1] ].y + this->Nodes[ Cell2Node[cell][2] ].y) / 3.0;

        this->CellVolume[cell] = 0.5 * fabs( this->Nodes[ Cell2Node[cell][0] ].x * ( this->Nodes[ Cell2Node[cell][1] ].y - this->Nodes[ Cell2Node[cell][2] ].y ) 
                                           + this->Nodes[ Cell2Node[cell][1] ].x * ( this->Nodes[ Cell2Node[cell][2] ].y - this->Nodes[ Cell2Node[cell][0] ].y ) 
                                           + this->Nodes[ Cell2Node[cell][2] ].x * ( this->Nodes[ Cell2Node[cell][0] ].y - this->Nodes[ Cell2Node[cell][1] ].y ) );
    }
}

void mshReader::computeFaceGeometry()
{
    this->FaceCenter.resize   ( Face2Node.size() );
    this->FaceNormal.resize   ( Face2Node.size() );
    this->FaceArea.resize     ( Face2Node.size() );
    this->FaceDistance.resize ( Face2Node.size() );

    for( int face = 0; face < this->Face2Node.size(); ++face )
    {
        this->FaceCenter[face].x  = 0.5 * ( this->Nodes[ this->Face2Node[face][0] ].x + this->Nodes[ this->Face2Node[face][1] ].x );
        this->FaceCenter[face].y  = 0.5 * ( this->Nodes[ this->Face2Node[face][0] ].y + this->Nodes[ this->Face2Node[face][1] ].y );

        this->FaceArea[face] = sqrt( ( this->Nodes[ this->Face2Node[face][1] ].x - this->Nodes[ this->Face2Node[face][0] ].x )
                                   * ( this->Nodes[ this->Face2Node[face][1] ].x - this->Nodes[ this->Face2Node[face][0] ].x ) 
                                   + ( this->Nodes[ this->Face2Node[face][1] ].y - this->Nodes[ this->Face2Node[face][0] ].y )
                                   * ( this->Nodes[ this->Face2Node[face][1] ].y - this->Nodes[ this->Face2Node[face][0] ].y ) );

        this->FaceNormal[face].x = - ( this->Nodes[ this->Face2Node[face][1] ].y - this->Nodes[ this->Face2Node[face][0] ].y ) / this->FaceArea[face];
        this->FaceNormal[face].y =   ( this->Nodes[ this->Face2Node[face][1] ].x - this->Nodes[ this->Face2Node[face][0] ].x ) / this->FaceArea[face];

        this->FaceDistance[face] = 0.0;
        if( this->Face2Cell[face][0] != -1 ) this->FaceDistance[face] += this->normalDistanceFace2Cell(face, Face2Cell[face][0]);
        if( this->Face2Cell[face][1] != -1 ) this->FaceDistance[face] += this->normalDistanceFace2Cell(face, Face2Cell[face][1]);
    }
}

void mshReader::computeCellMinDx()
{
    this->CellMinDx.resize( this->Cell2Node.size() );

    for( int cell = 0; cell < this->Cell2Node.size(); ++cell )
    {
        this->CellMinDx[cell] = 1.0e99;

        double distance;
        distance = fabs( this->FaceNormal[ Cell2Face[cell][0] ].x * this->Nodes[ Cell2Node[cell][2] ].x
                       + this->FaceNormal[ Cell2Face[cell][0] ].y * this->Nodes[ Cell2Node[cell][2] ].y );
        this->CellMinDx[cell] = min( this->CellMinDx[cell], distance );

        double distance;
        distance = fabs( this->FaceNormal[ Cell2Face[cell][1] ].x * this->Nodes[ Cell2Node[cell][0] ].x
                       + this->FaceNormal[ Cell2Face[cell][1] ].y * this->Nodes[ Cell2Node[cell][0] ].y );
        this->CellMinDx[cell] = min( this->CellMinDx[cell], distance );

        double distance;
        distance = fabs( this->FaceNormal[ Cell2Face[cell][2] ].x * this->Nodes[ Cell2Node[cell][1] ].x
                       + this->FaceNormal[ Cell2Face[cell][2] ].y * this->Nodes[ Cell2Node[cell][1] ].y );
        this->CellMinDx[cell] = min( this->CellMinDx[cell], distance );
    }
}

bool mshReader::findPeriodicCells()
{
    for(int left = 0; left < Face2Node.size(); ++left)
    {
        if( Face2BC[left] == -1 ) continue;
        if( BCs[ Face2BC[left] ] == periodicLeft && ( Face2Cell[left][0] != -1 || Face2Cell[left][1] != -1 ) )
        {
            for(int right = 0; right < Face2Node.size(); ++right)
            {
                if( BCs[ Face2BC[right] ] == periodicRight )
                {
                    Vec2 connection;

                    connection.x = this->FaceCenter[right].x - this->FaceCenter[left].x;
                    connection.y = this->FaceCenter[right].y - this->FaceCenter[left].y;

                    double projection = this->FaceNormal[left].x * connection.x + this->FaceNormal[left].y * connection.y;
                    double distance = sqrt( connection.x * connection.x + connection.y * connection.y );

                    if( fabs( projection - distance ) < 1.0e-12 )
                    {
                        int leftCell, leftEmpty;
                        int rightCell, rightEmpty;

                        if      (Face2Cell[left][0] != -1)
                        {
                            leftCell  = Face2Cell[left][0];
                            leftEmpty = 0;
                        }
                        else if (Face2Cell[left][1] != -1)
                        {
                            leftCell  = Face2Cell[left][1];
                            leftEmpty = 1;
                        }
                        else
                        {
                            cout << "Error: Left Cell already has two Neighbors" << endl; return false;
                        }

                        if      (Face2Cell[right][0] != -1)
                        {
                            rightCell  = Face2Cell[right][0];
                            rightEmpty = 0;
                        }
                        else if (Face2Cell[right][1] != -1)
                        {
                            rightCell  = Face2Cell[right][1];
                            rightEmpty = 1;
                        }
                        else
                        {
                            cout << "Error: Right Cell already has two Neighbors" << endl; return false;
                        }

                        Face2Cell[left ][leftEmpty ] = rightCell;
                        Face2Cell[right][rightEmpty] = leftCell;

                        double distance = this->FaceDistance[left] + this->FaceDistance[right];
                        this->FaceDistance[left ] = distance;
                        this->FaceDistance[right] = distance;

                        break;
                    }
                }
            }

            if( Face2Cell[left][0] != -1 || Face2Cell[left][0] != -1 )
            {
                cout << "Error: No periodic matching Cell found" << endl; return false;
            }
        }
    }

    return true;
}

double mshReader::normalDistanceFace2Cell(int face, int cell)
{
    return fabs( ( this->CellCenter[cell].x - this->FaceCenter[face].x ) * this->FaceNormal[face].x
               + ( this->CellCenter[cell].y - this->FaceCenter[face].y ) * this->FaceNormal[face].y );
}
