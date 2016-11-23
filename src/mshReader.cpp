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
    cout << "Start reading: " << filename << endl;
    ifstream file;
    file.open(filename);

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return false;
    }

    string buffer;

    getline(file, buffer);
    getline(file, buffer);
    getline(file, buffer);

    if( ! this->readPhysicalNames(file) ) return false;

    if( ! this->readNodes(file) )         return false;

    if( ! this->NodeCheck() )             return false;

    if( ! this->readElements(file) )      return false;
    
    this->computeCellGeometry();
    this->computeFaceGeometry();
    this->computeCellMinDx();

    if( ! this->findPeriodicCells() )     return false;

    this->createGhostCells();

    if( ! this->FaceCheck() )             return false;

    cout << "Complete Mesh read!" << endl;

    return true;
}

bool mshReader::readPhysicalNames(ifstream& file)
{
    cout << "Read Boundary Conditions:";

    string buffer;
    getline(file, buffer);
    if( buffer.compare("$PhysicalNames") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }

    idType nPhysicalNames;
    file >> nPhysicalNames;
    getline(file, buffer); // go to next line

    for( idType i = 0; i < nPhysicalNames; ++i )
    {
        getline(file, buffer);
        stringstream bufferStream(buffer);
        
        idType PhysicalNameType, ID;
        string type;

        bufferStream >> PhysicalNameType;
        if(PhysicalNameType != 1) continue;

        bufferStream >> ID >> type;

        this->BCIDs.push_back(ID);

        type = type.substr(1, type.length()-2 );

        if     ( type.compare("periodic") == 0 ) this->BCs.push_back(periodic);
        else if( type.compare("wall")     == 0 ) this->BCs.push_back(wall);
        else if( type.compare("inlet")    == 0 ) this->BCs.push_back(inlet);
        else if( type.compare("outlet")   == 0 ) this->BCs.push_back(outlet);
        else { cout << "Error: Wrong BC" << endl; return false; }
    }

    getline(file, buffer);
    if( buffer.compare("$EndPhysicalNames") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }
    
    cout << " done!" << endl;

    return true;
}

bool mshReader::readNodes(ifstream & file)
{
    cout << "Read Node Coordinates:";

    string buffer;
    getline(file, buffer);
    if( buffer.compare("$Nodes") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }

    idType nNodes;
    file >> nNodes;
    getline(file, buffer); // go to next line

    for( idType i = 0; i < nNodes; ++i )
    {
        getline(file, buffer);
        stringstream bufferStream(buffer);
        
        idType ID;
        Vec2 tmpNode(0,0);

        bufferStream >> ID;
        if(ID != i+1) { cout << "Error: Invalid Node ID" << endl; return false; }

        bufferStream >> tmpNode.x >> tmpNode.y;

        //cout << endl << "Read Node " << this->Nodes.size() << " at ( " << tmpNode.x << ", " << tmpNode.y << " )" << endl;

        this->Nodes.push_back(tmpNode);
    }

    getline(file, buffer);
    if( buffer.compare("$EndNodes") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }
    
    cout << " done!" << endl;

    return true;
}

bool mshReader::readElements(ifstream & file)
{
    cout << "Read Cells and Boundary Interfaces:";

    string buffer;
    getline(file, buffer);
    if( buffer.compare("$Elements") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }

    idType nElements;
    file >> nElements;
    getline(file, buffer); // go to next line

    for( idType i = 0; i < nElements; ++i )
    {
        getline(file, buffer);
        stringstream bufferStream(buffer);
        
        idType ID, type;

        bufferStream >> ID >> type;

        if (type == 1) if( ! newFace(buffer)       ) return false;
        if (type == 2) if( ! newCell(buffer, tri ) ) return false;
        if (type == 3) if( ! newCell(buffer, quad) ) return false;
    }

    getline(file, buffer);
    if( buffer.compare("$EndElements") != 0 ) { cout << "Error: Wrong Key" << endl; return false; }
    
    cout << " done!" << endl;

    return true;
}

bool mshReader::newFace(string buffer)
{
    stringstream bufferStream(buffer);

    idType ID;
    idType BCID;
    idType N1, N2;
    idType tmp;

    bufferStream >> ID >> tmp >> tmp >> BCID >> tmp >> N1 >> N2;
    
    for(idType i = 0; i < this->Face2Node.size(); ++i)
        if (  ( this->Face2Node[i][0] == N1-1 && this->Face2Node[i][1] == N2-1 )
           || ( this->Face2Node[i][1] == N1-1 && this->Face2Node[i][0] == N2-1 ) )
        {
            this->Face2BC[i] = findIndex( this->BCIDs, BCID );
            return true;
        }

    array<idType,2> tmpFace2Node;
    tmpFace2Node[0] = N1-1;
    tmpFace2Node[1] = N2-1;
    this->Face2Node.push_back(tmpFace2Node);

    array<idType,2> tmpFace2Cell;
    tmpFace2Cell[0] = -1;
    tmpFace2Cell[1] = -1;
    this->Face2Cell.push_back(tmpFace2Cell);

    array<bool,2> tmpFace2CellAdd = {false, false};
    this->Face2CellAdd.push_back(tmpFace2CellAdd);
    
    this->Face2BC.push_back( findIndex( this->BCIDs, BCID ) );

    return true;
}

bool mshReader::newCell(string buffer, CellType type)
{
    stringstream bufferStream(buffer);

    idType ID;
    idType BCID;
    idType N1, N2, N3, N4;
    idType tmp;
    array<idType,4> tmpCell2Node;

    if( type == tri )
    {
        bufferStream >> ID >> tmp >> tmp >> tmp >> tmp >> N1 >> N2 >> N3;

        tmpCell2Node[0] = N1-1;
        tmpCell2Node[1] = N2-1;
        tmpCell2Node[2] = N3-1;
        tmpCell2Node[3] = -1;
    }
    else if (type == quad)
    {
        bufferStream >> ID >> tmp >> tmp >> tmp >> tmp >> N1 >> N2 >> N3 >> N4;

        tmpCell2Node[0] = N1-1;
        tmpCell2Node[1] = N2-1;
        tmpCell2Node[2] = N3-1;
        tmpCell2Node[3] = N4-1;
    }

    array<idType, 4> tmpCell2Face;
    tmpCell2Face[0] = -1;
    tmpCell2Face[1] = -1;
    tmpCell2Face[2] = -1;
    tmpCell2Face[3] = -1;

    this->linkExistingFaces( tmpCell2Node, tmpCell2Face, this->Cell2Node.size(), type);

    this->createMissingFaces( tmpCell2Node, tmpCell2Face, this->Cell2Node.size(), type);

    this->Cell2Type.push_back(type);
    this->Cell2Node.push_back(tmpCell2Node);
    this->Cell2Face.push_back(tmpCell2Face);

    this->Cell2BC.push_back(-1);

    return true;
}

bool mshReader::linkExistingFaces(array<idType, 4> tmpCell2Node, array<idType, 4>& tmpCell2Face, idType CellID, CellType type)
{
    idType nNodes;
    if     (type == tri)  nNodes = 3;
    else if(type == quad) nNodes = 4;

    for( idType face = 0; face < Face2Node.size(); ++face )
    {
        for( idType cellNode = 0; cellNode < nNodes; ++cellNode )
        {
            if      ( this->Face2Node[face][0] == tmpCell2Node[ cellNode%nNodes ] && this->Face2Node[face][1] == tmpCell2Node[ (cellNode+1)%nNodes ] )
            {
                tmpCell2Face[cellNode] = face;
                this->Face2Cell   [face][1] = CellID;
                this->Face2CellAdd[face][1] = true;
            }
            else if ( this->Face2Node[face][1] == tmpCell2Node[ cellNode%nNodes ] && this->Face2Node[face][0] == tmpCell2Node[ (cellNode+1)%nNodes ] )
            {
                tmpCell2Face[cellNode] = face;
                this->Face2Cell   [face][0] = CellID;
                this->Face2CellAdd[face][0] = true;
            }  
        }
    }

    return true;
}

bool mshReader::createMissingFaces(array<idType, 4> tmpCell2Node, array<idType, 4>& tmpCell2Face, idType CellID, CellType type)
{
    idType nFaces;
    if     (type == tri)   nFaces = 3;
    else if(type == quad) nFaces = 4;

    for( idType cellFace = 0; cellFace < nFaces; ++cellFace )
    {
        if(tmpCell2Face[cellFace] == -1)
        {
            array<idType,2> tmpFace2Node;
            tmpFace2Node[0] = tmpCell2Node[ cellFace   %nFaces];
            tmpFace2Node[1] = tmpCell2Node[(cellFace+1)%nFaces];
            Face2Node.push_back(tmpFace2Node);

            array<idType,2> tmpFace2Cell;
            tmpFace2Cell[0] = -1;
            tmpFace2Cell[1] = CellID;
            Face2Cell.push_back(tmpFace2Cell);

            array<bool,2> tmpFace2CellAdd = {true, true};
            this->Face2CellAdd.push_back(tmpFace2CellAdd);

            Face2BC.push_back(-1);

            tmpCell2Face[cellFace] = Face2Cell.size()-1;
        }  
    }
    
    return true;
}

void mshReader::computeCellGeometry()
{
    cout << "Compute Cell geometry:";
    this->CellCenter.resize( this->Cell2Node.size() );
    this->CellVolume.resize( this->Cell2Node.size() );

    for( idType cell = 0; cell < this->Cell2Node.size(); ++cell )
    {
        if      ( this->Cell2Type[cell] == tri )
        {
            this->CellCenter[cell].x =  (this->Nodes[ Cell2Node[cell][0] ].x + this->Nodes[ Cell2Node[cell][1] ].x + this->Nodes[ Cell2Node[cell][2] ].x) / 3.0;
            this->CellCenter[cell].y =  (this->Nodes[ Cell2Node[cell][0] ].y + this->Nodes[ Cell2Node[cell][1] ].y + this->Nodes[ Cell2Node[cell][2] ].y) / 3.0;

            this->CellVolume[cell] = 0.5 * fabs( this->Nodes[ Cell2Node[cell][0] ].x * ( this->Nodes[ Cell2Node[cell][1] ].y - this->Nodes[ Cell2Node[cell][2] ].y ) 
                                               + this->Nodes[ Cell2Node[cell][1] ].x * ( this->Nodes[ Cell2Node[cell][2] ].y - this->Nodes[ Cell2Node[cell][0] ].y ) 
                                               + this->Nodes[ Cell2Node[cell][2] ].x * ( this->Nodes[ Cell2Node[cell][0] ].y - this->Nodes[ Cell2Node[cell][1] ].y ) );
        }
        else if ( this->Cell2Type[cell] == quad )
        {
            Vec2 triCenter[2];
            triCenter[0].x =  (this->Nodes[ Cell2Node[cell][0] ].x + this->Nodes[ Cell2Node[cell][1] ].x +                                       this->Nodes[ Cell2Node[cell][3] ].x) / 3.0;
            triCenter[0].y =  (this->Nodes[ Cell2Node[cell][0] ].y + this->Nodes[ Cell2Node[cell][1] ].y +                                       this->Nodes[ Cell2Node[cell][3] ].y) / 3.0;
            triCenter[1].y =  (                                      this->Nodes[ Cell2Node[cell][1] ].y + this->Nodes[ Cell2Node[cell][2] ].y + this->Nodes[ Cell2Node[cell][3] ].y) / 3.0;
            triCenter[1].x =  (                                      this->Nodes[ Cell2Node[cell][1] ].x + this->Nodes[ Cell2Node[cell][2] ].x + this->Nodes[ Cell2Node[cell][3] ].x) / 3.0;

            double triVolume[2];
            triVolume[0] = 0.5 * fabs( this->Nodes[ Cell2Node[cell][0] ].x * ( this->Nodes[ Cell2Node[cell][1] ].y - this->Nodes[ Cell2Node[cell][3] ].y ) 
                                     + this->Nodes[ Cell2Node[cell][1] ].x * ( this->Nodes[ Cell2Node[cell][3] ].y - this->Nodes[ Cell2Node[cell][0] ].y ) 
                                     + this->Nodes[ Cell2Node[cell][3] ].x * ( this->Nodes[ Cell2Node[cell][0] ].y - this->Nodes[ Cell2Node[cell][1] ].y ) );
            triVolume[1] = 0.5 * fabs( this->Nodes[ Cell2Node[cell][2] ].x * ( this->Nodes[ Cell2Node[cell][3] ].y - this->Nodes[ Cell2Node[cell][1] ].y ) 
                                     + this->Nodes[ Cell2Node[cell][3] ].x * ( this->Nodes[ Cell2Node[cell][1] ].y - this->Nodes[ Cell2Node[cell][2] ].y ) 
                                     + this->Nodes[ Cell2Node[cell][1] ].x * ( this->Nodes[ Cell2Node[cell][2] ].y - this->Nodes[ Cell2Node[cell][3] ].y ) );

            this->CellVolume[cell]   = triVolume[0] + triVolume[1];
            this->CellCenter[cell].x = ( triCenter[0].x * triVolume[0] + triCenter[1].x * triVolume[1] ) / this->CellVolume[cell];
            this->CellCenter[cell].y = ( triCenter[0].y * triVolume[0] + triCenter[1].y * triVolume[1] ) / this->CellVolume[cell];
        }
    }
    cout << " done!" << endl;
}

void mshReader::computeCellLSCoeff()
{
    cout << "Compute Cell Least Square Coefficients:";
    this->CellLSCoeff.resize( this->Cell2Node.size() );

    for( idType cell = 0; cell < this->Cell2Node.size(); ++cell )
    {
        idType nNeighbors;
        if      ( this->Cell2Type[cell] == tri  ) nNeighbors = 3;
        else if ( this->Cell2Type[cell] == quad ) nNeighbors = 4;

        for( idType face = 0; face < nNeighbors; ++face )
        {
            double dx = this->CellCenter[ this->getNeighborCell( this->Cell2Face[cell][face], cell) ].x - this->CellCenter[cell].x;
            double dy = this->CellCenter[ this->getNeighborCell( this->Cell2Face[cell][face], cell) ].y - this->CellCenter[cell].y;

            double distance = sqrt(dx*dx + dy*dy);

            this->CellLSCoeff[cell][0] += dx*dx;
            this->CellLSCoeff[cell][1] += dx*dy;
            this->CellLSCoeff[cell][2] += dy*dy;
        }

        this->CellLSCoeff[cell][0] = sqrt( this->CellLSCoeff[cell][0] );
        this->CellLSCoeff[cell][1] = this->CellLSCoeff[cell][1] / this->CellLSCoeff[cell][0];
        this->CellLSCoeff[cell][2] = sqrt( this->CellLSCoeff[cell][2] - this->CellLSCoeff[cell][1]*this->CellLSCoeff[cell][1] );
    }

    cout << " done!" << endl;
}

void mshReader::computeFaceGeometry()
{
    cout << "Compute Face gemometry:";
    this->FaceCenter.resize   ( Face2Node.size() );
    this->FaceNormal.resize   ( Face2Node.size() );
    this->FaceArea.resize     ( Face2Node.size() );
    this->FaceDistance.resize ( Face2Node.size() );

    for( idType face = 0; face < this->Face2Node.size(); ++face )
    {
        this->FaceCenter[face].x  = 0.5 * ( this->Nodes[ this->Face2Node[face][0] ].x + this->Nodes[ this->Face2Node[face][1] ].x );
        this->FaceCenter[face].y  = 0.5 * ( this->Nodes[ this->Face2Node[face][0] ].y + this->Nodes[ this->Face2Node[face][1] ].y );

        this->FaceArea[face] = sqrt( ( this->Nodes[ this->Face2Node[face][1] ].x - this->Nodes[ this->Face2Node[face][0] ].x )
                                   * ( this->Nodes[ this->Face2Node[face][1] ].x - this->Nodes[ this->Face2Node[face][0] ].x ) 
                                   + ( this->Nodes[ this->Face2Node[face][1] ].y - this->Nodes[ this->Face2Node[face][0] ].y )
                                   * ( this->Nodes[ this->Face2Node[face][1] ].y - this->Nodes[ this->Face2Node[face][0] ].y ) );
        
        //                  Compute Normal
        // ========================================================================
        //      -----[1]-------
        //            |             The normal is computed such that it points
        //        n   |             to the right when the one follows the
        //      <-----|             vector  from the first to the second node.
        //    posCell | negCell     
        //            |             n =  (0 0 1) x (N1 - N0)
        //      -----[0]-------    
        // ========================================================================
        this->FaceNormal[face].x = - ( this->Nodes[ this->Face2Node[face][1] ].y - this->Nodes[ this->Face2Node[face][0] ].y ) / this->FaceArea[face];
        this->FaceNormal[face].y =   ( this->Nodes[ this->Face2Node[face][1] ].x - this->Nodes[ this->Face2Node[face][0] ].x ) / this->FaceArea[face];

        this->FaceDistance[face] = 0.0;
        if( this->Face2Cell[face][0] != -1 ) this->FaceDistance[face] += this->normalDistanceFace2Cell(face, Face2Cell[face][0]);
        if( this->Face2Cell[face][1] != -1 ) this->FaceDistance[face] += this->normalDistanceFace2Cell(face, Face2Cell[face][1]);
    }
    cout << " done!" << endl;
}

void mshReader::computeCellMinDx()
{
    cout << "Compute minimal Cell extend:";
    this->CellMinDx.resize( this->Cell2Node.size() );

    for( idType cell = 0; cell < this->Cell2Node.size(); ++cell )
    {
        this->CellMinDx[cell] = 1.0e99;

        idType nNodes;
        if      ( this->Cell2Type[cell] == tri  ) nNodes = 3;
        else if ( this->Cell2Type[cell] == quad ) nNodes = 4;

        for(idType cellFace = 0; cellFace < nNodes; ++cellFace)          // loop over faces
        {
            for( idType cellNode = 0; cellNode < nNodes; ++cellNode)     // loop over nodes not on the interface
            {
                if(  this->Cell2Node[cell][cellNode] != this->Face2Node[ this->Cell2Face[cell][cellFace] ][0] 
                  && this->Cell2Node[cell][cellNode] != this->Face2Node[ this->Cell2Face[cell][cellFace] ][1] )
                {
                    this->CellMinDx[cell] = min( this->CellMinDx[cell], this->normalDistanceFace2Point( this->Cell2Face[cell][cellFace], this->Nodes[ this->Cell2Node[cell][cellNode] ] ) );
                }
            }
        }

    }
    cout << " done!" << endl;
}

bool mshReader::findPeriodicCells()
{
    cout << "Find periodic Cell:";
    for(idType left = 0; left < Face2Node.size(); ++left)
    {
        if( this->Face2BC[left] == -1 ) continue;
        if( this->BCs[ Face2BC[left] ] == periodic && ( this->Face2Cell[left][0] != -1 || this->Face2Cell[left][1] != -1 ) )
        {
            for(idType right = 0; right < Face2Node.size(); ++right)
            {
                if( BCs[ Face2BC[right] ] == periodic && left != right )
                {
                    Vec2 connection;

                    connection.x = this->FaceCenter[right].x - this->FaceCenter[left].x;
                    connection.y = this->FaceCenter[right].y - this->FaceCenter[left].y;

                    double projection = fabs( this->FaceNormal[left].x * connection.x + this->FaceNormal[left].y * connection.y );
                    double distance = sqrt( connection.x * connection.x + connection.y * connection.y );

                    if( fabs( projection - distance ) < 1.0e-10 )
                    {
                        idType leftCell, leftEmpty;
                        idType rightCell, rightEmpty;

                        if      (this->Face2Cell[left][0] != -1)
                        {
                            leftCell  = this->Face2Cell[left][0];
                            leftEmpty = 1;
                        }
                        else if (this->Face2Cell[left][1] != -1)
                        {
                            leftCell  = this->Face2Cell[left][1];
                            leftEmpty = 0;
                        }
                        else
                        {
                            cout << "Error: Left Cell already has two Neighbors" << endl; return false;
                        }

                        if      (this->Face2Cell[right][0] != -1)
                        {
                            rightCell  = this->Face2Cell[right][0];
                            rightEmpty = 1;
                        }
                        else if (this->Face2Cell[right][1] != -1)
                        {
                            rightCell  = this->Face2Cell[right][1];
                            rightEmpty = 0;
                        }
                        else
                        {
                            cout << "Error: Right Cell already has two Neighbors" << endl; return false;
                        }

                        this->Face2Cell[left ][leftEmpty ] = rightCell;
                        this->Face2Cell[right][rightEmpty] = leftCell;

                        double distance = this->FaceDistance[left] + this->FaceDistance[right];
                        this->FaceDistance[left ] = distance;
                        this->FaceDistance[right] = distance;

                        break;
                    }
                }
            }

            if( Face2Cell[left][0] == -1 || Face2Cell[left][0] == -1 )
            {
                cout << "Error: No periodic matching Cell found for face " << left << endl; return false;
            }
        }
    }

    cout << " done!" << endl;
    return true;
}

void mshReader::createGhostCells()
{
    for( idType face = 0; face < this->Face2Node.size(); ++face )
    {
        idType emptyFace2Cell;
        if      (this->Face2Cell[face][0] == -1) emptyFace2Cell = 0;
        else if (this->Face2Cell[face][1] == -1) emptyFace2Cell = 1;
        else continue;

        Vec2 vector( ( (emptyFace2Cell == 1)?(3.0):(-3.0) ) * this->FaceDistance[face] * this->FaceNormal[face].x,
                     ( (emptyFace2Cell == 1)?(3.0):(-3.0) ) * this->FaceDistance[face] * this->FaceNormal[face].y );

        this->Nodes.push_back( Vec2( this->FaceCenter[face].x + vector.x,
                                     this->FaceCenter[face].y + vector.y ) );

        array<idType,4> newCell2Node;
        if( emptyFace2Cell == 0 )
        {
            newCell2Node[0] = this->Face2Node[face][0];
            newCell2Node[1] = this->Face2Node[face][1];
        }
        else
        {
            newCell2Node[0] = this->Face2Node[face][1];
            newCell2Node[1] = this->Face2Node[face][0];
        }
        newCell2Node[2] = this->Nodes.size()-1;
        newCell2Node[3] = -1;

        array<idType,4> newCell2Face;
        newCell2Face[0] = face;
        newCell2Face[1] = -1;
        newCell2Face[2] = -1;
        newCell2Face[3] = -1;

        this->Cell2Type.push_back(tri);
        this->Cell2Node.push_back(newCell2Node);
        this->Cell2Face.push_back(newCell2Face);
        this->Cell2BC.push_back(this->Face2BC[face]);

        this->CellCenter.push_back(Vec2(0.0,0.0));
        this->CellVolume.push_back(0.0);
        this->CellMinDx.push_back(0.0);

        this->Face2Cell[face][emptyFace2Cell] = Cell2Node.size()-1;
        this->FaceDistance[face] *= 2.0;
    }
}

bool mshReader::NodeCheck()
{
    for( idType node_1 = 0; node_1 < this->Nodes.size(); ++node_1 )
    for( idType node_2 = 0; node_2 < this->Nodes.size(); ++node_2 )
    {
        if  (  node_1 != node_2 
            && fabs( Nodes[node_1].x - Nodes[node_2].x ) < 1.0e-10
            && fabs( Nodes[node_1].y - Nodes[node_2].y ) < 1.0e-10 )
        {
            cout << "Error: Identical Nodes found: " << node_1 << " and " << node_2 << endl;
            cout << "dX = " << Nodes[node_1].x << ", " << Nodes[node_2].x << ", Dy = " << Nodes[node_1].y << ", " << Nodes[node_2].y << endl;
            return false;
        }
    }
    return true;
}

bool mshReader::FaceCheck()
{
    for( idType face = 0; face < this->Face2Node.size(); ++face )
    {
        if( this->Face2Cell[face][0] == -1 || this->Face2Cell[face][1] == -1 )
        {
            cout << "Error: Free cell faces detected at face " << face << endl;
            return false;
        }
    }
    return true;
}

double mshReader::normalDistanceFace2Cell(idType face, idType cell)
{
    return fabs( ( this->CellCenter[cell].x - this->FaceCenter[face].x ) * this->FaceNormal[face].x
               + ( this->CellCenter[cell].y - this->FaceCenter[face].y ) * this->FaceNormal[face].y );
}

double mshReader::normalDistanceFace2Point(idType face, Vec2 point)
{
    return fabs( ( point.x - this->FaceCenter[face].x ) * this->FaceNormal[face].x
               + ( point.y - this->FaceCenter[face].y ) * this->FaceNormal[face].y );
}

idType mshReader::getNeighborCell(idType face, idType askingCell)
{
    if      ( this->Face2Cell[face][0] == askingCell ) return this->Face2Cell[face][1];
    else if ( this->Face2Cell[face][1] == askingCell ) return this->Face2Cell[face][0];
    
    return -1;
}


