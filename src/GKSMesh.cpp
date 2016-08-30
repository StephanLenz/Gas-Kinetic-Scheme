
#include "GKSMesh.h"
#include "Cell.h"
#include "CompressibleInterface.h"
#include "Types.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>    //min()
#include <chrono>

# define M_PI           3.14159265358979323846  /* pi */

using namespace std;

GKSMesh::GKSMesh()
{
}

GKSMesh::GKSMesh(Parameters param, FluidParameter fluidParam)
{
    this->param = param;
    this->fluidParam = fluidParam;
    this->iter = 0;
}


GKSMesh::~GKSMesh()
{
}

void GKSMesh::generateRectMeshGraded(InterfaceType type, double lengthX, double lengthY, int nx, int ny, double gradingX, double gradingY)
{

    this->lengthX = lengthX;
    this->lengthY = lengthY;

    // make nx and ny an even number
    nx = 2 * (nx/2);
    ny = 2 * (ny/2);
    
    double etaX = pow( gradingX, 1.0 / (nx/2 - 1) );
    double etaY = pow( gradingY, 1.0 / (ny/2 - 1) );

    double dx0;
    if( fabs(etaX - 1.0) > 1.0e-12 )
        dx0 = 0.5*this->lengthX * (1-etaX)/(1-pow(etaX, nx/2));
    else
        dx0 = this->lengthX / double(nx);

    double dy0;
    if( fabs(etaY - 1.0) > 1.0e-12 )
        dy0 = 0.5*this->lengthY * (1-etaY)/(1-pow(etaY, ny/2));
    else
        dy0 = this->lengthY / double(ny);

    //=========================================================================
    //=========================================================================
    //		Computation of the coordinates and spacings
    //=========================================================================
    //=========================================================================

    double* CellSpacingsX = new double[nx];
    for(int i = 0; i < nx/2; i++)
    {
        CellSpacingsX[nx/2 - i - 1] = dx0 * pow( etaX, i);
        CellSpacingsX[nx/2 + i    ] = dx0 * pow( etaX, i);
    }

    double* CellSpacingsY = new double[ny];
    for(int i = 0; i < ny/2; i++)
    {
        CellSpacingsY[ny/2 - i - 1] = dy0 * pow( etaY, i);
        CellSpacingsY[ny/2 + i    ] = dy0 * pow( etaY, i);
    }

    // ========================================================================
    
    double* NodesX = new double[nx+1]; 
    double* NodesY = new double[ny+1];

    double sumX = 0.0;
    for(int i = 0; i < nx+1; i++){
        NodesX[i] = sumX;
        sumX += CellSpacingsX[i];
    }
   
    double sumY = 0.0;
    for(int i = 0; i < ny+1; i++){
        NodesY[i] = sumY;
        sumY += CellSpacingsY[i];
    }
    
	//=========================================================================
	//=========================================================================
	//		Node generation
	//=========================================================================
	//=========================================================================

    double heightDiff = 0.5;

    for (int i = 0; i < ny + 1; i++)       // Y-Direction
	{
		for (int j = 0; j < nx + 1; j++)   // X-Direction
		{
            //float2* tmpNode = new float2( NodesX[j], NodesY[i] );
            //float2* tmpNode = new float2( NodesX[j], NodesY[i] + NodesX[j] / this->lengthX * heightDiff );

            float2* tmpNode = new float2( NodesX[j], NodesY[i] - 1.0 * (NodesX[j] - this->lengthX)*NodesX[j]*(NodesY[i] - this->lengthY)*NodesY[i] );

            //float2* tmpNode = new float2();
            //tmpNode->x = NodesX[j];
            ////tmpNode->y = NodesY[i] - 0.25 * (NodesX[j] - this->lengthX)*NodesX[j] * sin( (NodesY[i] - 0.5*this->lengthY) * 2.0 * M_PI/this->lengthY );
            //tmpNode->y = NodesY[i] - 0.05 * sin( NodesX[j] * 2.0 * M_PI/this->lengthX ) * sin( (NodesY[i] - 0.5*this->lengthY) * 2.0 * M_PI/this->lengthY );

			this->NodeList.push_back(tmpNode);
		}
	}
    
    
	//=========================================================================
	//=========================================================================
	//		Node Rotation
	//=========================================================================
	//=========================================================================

    double angle = 0.0;//0.25*M_PI;

    for (vector<float2*>::iterator i = NodeList.begin(); i != NodeList.end(); ++i)
	{
        double x = (*i)->x;
        double y = (*i)->y;
        (*i)->x = cos(angle) * x - sin(angle) * y;
        (*i)->y = sin(angle) * x + cos(angle) * y;
	}

	//=========================================================================
	//=========================================================================
	//		Cell generation
	//=========================================================================
	//=========================================================================
	for (int i = 0; i < ny; i++)       // Y-Direction
	{
		for (int j = 0; j < nx; j++)   // X-Direction
		{
            float2* tmpNodes[4];
            tmpNodes[0] = this->NodeList[(i+0)*(nx+1) + j + 0]; // bottom left
            tmpNodes[1] = this->NodeList[(i+0)*(nx+1) + j + 1]; // bottom right
            tmpNodes[2] = this->NodeList[(i+1)*(nx+1) + j + 1]; // top right
            tmpNodes[3] = this->NodeList[(i+1)*(nx+1) + j + 0]; // top left

			Cell* tmpCell = new Cell(type, tmpNodes, NULL, this->fluidParam);
			// add interface to list
			this->CellList.push_back(tmpCell);
		}
	}

	//=========================================================================
	//=========================================================================
	//						F interface generation
	//=========================================================================
	//=========================================================================
	for (int i = 0; i < ny; i++)       // Y-Direction
	{
		for (int j = 0; j < nx + 1; j++)    // X-Direction
		{
            float2* tmpNodes[2];
            tmpNodes[0] = this->NodeList[(i+1)*(nx+1) + j]; // top
            tmpNodes[1] = this->NodeList[(i+0)*(nx+1) + j]; // bottom

            Cell* negCell = NULL;
            Cell* posCell = NULL;

            if( j != 0  ) 
                negCell = this->CellList[i*nx + (j  - 1)];
            else if( this->BoundaryConditionList[0]->getType() == periodic )
                negCell = this->CellList[i*nx + (nx - 1)];

            if( j != nx ) 
                posCell = this->CellList[i*nx + (j  - 0)];
            else if( this->BoundaryConditionList[2]->getType() == periodic )       
                posCell = this->CellList[i*nx + (0  - 0)];

            BoundaryCondition* currentBC = NULL;
            if( j == 0  && this->BoundaryConditionList[0]->getType() != periodic )
                currentBC = this->BoundaryConditionList[0];

            if( j == nx && this->BoundaryConditionList[2]->getType() != periodic )
                currentBC = this->BoundaryConditionList[2];
            
            bool posAdd = (j != nx);
            bool negAdd = (j != 0) ;

			// create a new interface with the adjacent cells
			Interface* tmpInterface = Interface::createInterface(type,negCell, posCell, negAdd, posAdd, tmpNodes, this->fluidParam, currentBC, lengthX);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}

	//=========================================================================
	//=========================================================================
	//						G interface generation
	//=========================================================================
	//=========================================================================
	for (int i = 0; i < ny + 1; i++)        // Y-Direction
	{
		for (int j = 0; j < nx; j++)   // X-Direction
		{
            float2* tmpNodes[2];
            tmpNodes[0] = this->NodeList[i*(nx+1) + (j+0)]; // left
            tmpNodes[1] = this->NodeList[i*(nx+1) + (j+1)]; // right

            Cell* negCell = NULL;
            Cell* posCell = NULL;

            if( i != 0  ) 
                negCell = this->CellList[(i -1)*nx + j];
            else if( this->BoundaryConditionList[1]->getType() == periodic )
                negCell = this->CellList[(ny-1)*nx + j];

            if( i != ny ) 
                posCell = this->CellList[(i -0)*nx + j];
            else if( this->BoundaryConditionList[3]->getType() == periodic )       
                posCell = this->CellList[(0 -0)*nx + j];

            BoundaryCondition* currentBC = NULL;
            if( i == 0  && this->BoundaryConditionList[1]->getType() != periodic ) 
                currentBC = this->BoundaryConditionList[1];

            if( i == ny && this->BoundaryConditionList[3]->getType() != periodic ) 
                currentBC = this->BoundaryConditionList[3];
            
            bool posAdd = (i != ny);
            bool negAdd = (i != 0) ;

			// create a new interface with the adjacent cells
			Interface* tmpInterface = Interface::createInterface(type,negCell, posCell, negAdd, posAdd, tmpNodes, this->fluidParam, currentBC, lengthY);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}

	//=========================================================================
	//=========================================================================
	//						create Ghostcells
	//=========================================================================
	//=========================================================================
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        if( (*i)->isBoundaryInterface() )
        {
            float2* tmpNodes[4];
            float2  normal = (*i)->getScaledNormal();
            
            tmpNodes[1] = (*i)->getNode(0);
            tmpNodes[2] = (*i)->getNode(1);

            tmpNodes[0] = new float2( tmpNodes[1]->x + 2.0*normal.x, tmpNodes[1]->y + 2.0*normal.y );
            tmpNodes[3] = new float2( tmpNodes[2]->x + 2.0*normal.x, tmpNodes[2]->y + 2.0*normal.y );

            this->NodeList.push_back(tmpNodes[0]);
            this->NodeList.push_back(tmpNodes[3]);

            Cell* tmpCell = new Cell(type, tmpNodes, (*i)->getBoundaryCondition(), this->fluidParam);

            tmpCell->addInterface(*i);
            (*i)->addCell(tmpCell);

            this->CellList.push_back(tmpCell);
        }
    }

	//=========================================================================
	//=========================================================================
	//						Compute MinDx
	//=========================================================================
	//=========================================================================
    for(vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        (*i)->computeMinDx();
    }
    
	//=========================================================================
    delete [] CellSpacingsX;
    delete [] NodesX;

    delete [] CellSpacingsY;
    delete [] NodesY;

	return;
}

void GKSMesh::generateRectMeshPeriodicGraded(InterfaceType type, double lengthX, double lengthY, int nx, int ny, double grading)
{

    this->lengthX = lengthX;
    this->lengthY = lengthY;

    // make ny an even number
    ny = 2 * (ny/2);
    
    double etaY = pow( grading, 1.0 / (ny/2 - 1) );

    double dx0 = this->lengthX / nx;

    double dy0;
    if( fabs(etaY - 1.0) > 1.0e-12 )
        dy0 = 0.5*this->lengthY * (1-etaY)/(1-pow(etaY, ny/2));
    else
        dy0 = this->lengthY / double(ny);

    //=========================================================================
    //=========================================================================
    //		Computation of the coordinates and spacings
    //=========================================================================
    //=========================================================================

    double* CellSpacingsX = new double[nx];
    double* CellSpacingsY = new double[ny];

    for(int i = 0; i < nx; i++){
        CellSpacingsX[i] = dx0;
    }

    for(int i = 0; i < ny/2; i++){
        CellSpacingsY[ny/2 - i - 1] = dy0 * pow( etaY, i);
        CellSpacingsY[ny/2 + i    ] = dy0 * pow( etaY, i);
    }

    // ========================================================================
    
    double* NodesX = new double[nx+1]; 
    double* NodesY = new double[ny+1];

    double sumX = 0.0;
    for(int i = 0; i < nx+1; i++){
        NodesX[i] = sumX;
        sumX += CellSpacingsX[i];
    }
   
    double sumY = 0.0;
    for(int i = 0; i < ny+1; i++){
        NodesY[i] = sumY;
        sumY += CellSpacingsY[i];
    }
    
	//=========================================================================
	//=========================================================================
	//		Node generation
	//=========================================================================
	//=========================================================================

    double heightDiff = 0.5;

    for (int i = 0; i < ny + 1; i++)       // Y-Direction
	{
		for (int j = 0; j < nx + 1; j++)   // X-Direction
		{
            //float2* tmpNode = new float2( NodesX[j], NodesY[i] );
            //float2* tmpNode = new float2( NodesX[j], NodesY[i] + NodesX[j] / this->lengthX * heightDiff );
            
            float2* tmpNode = new float2();
            tmpNode->x = NodesX[j];
            //tmpNode->y = NodesY[i] - 0.25 * (NodesX[j] - this->lengthX)*NodesX[j] * sin( (NodesY[i] - 0.5*this->lengthY) * 2.0 * M_PI/this->lengthY );
            tmpNode->y = NodesY[i] - 0.05 * sin( NodesX[j] * 2.0 * M_PI/this->lengthX ) * sin( (NodesY[i] - 0.5*this->lengthY) * 2.0 * M_PI/this->lengthY );

			this->NodeList.push_back(tmpNode);
		}
	}
    
    
	//=========================================================================
	//=========================================================================
	//		Node Rotation
	//=========================================================================
	//=========================================================================

    double angle = 0.0;//0.25*M_PI;

    for (vector<float2*>::iterator i = NodeList.begin(); i != NodeList.end(); ++i)       // Y-Direction
	{
        double x = (*i)->x;
        double y = (*i)->y;
        (*i)->x = cos(angle) * x - sin(angle) * y;
        (*i)->y = sin(angle) * x + cos(angle) * y;
	}

	//=========================================================================
	//=========================================================================
	//		Cell generation
	//=========================================================================
	//=========================================================================
	for (int i = 0; i < ny; i++)       // Y-Direction
	{
		for (int j = 0; j < nx; j++)   // X-Direction
		{
            float2* tmpNodes[4];
            tmpNodes[0] = this->NodeList[(i+0)*(nx+1) + j + 0]; // bottom left
            tmpNodes[1] = this->NodeList[(i+0)*(nx+1) + j + 1]; // bottom right
            tmpNodes[2] = this->NodeList[(i+1)*(nx+1) + j + 1]; // top right
            tmpNodes[3] = this->NodeList[(i+1)*(nx+1) + j + 0]; // top left

			Cell* tmpCell = new Cell(type, tmpNodes, NULL, this->fluidParam);
			// add interface to list
			this->CellList.push_back(tmpCell);
		}
	}

	//=========================================================================
	//=========================================================================
	//						F interface generation
	//=========================================================================
	//=========================================================================
	for (int i = 0; i < ny; i++)       // Y-Direction
	{
		for (int j = 0; j < nx + 1; j++)    // X-Direction
		{
            float2* tmpNodes[2];
            tmpNodes[0] = this->NodeList[(i+1)*(nx+1) + j]; // top
            tmpNodes[1] = this->NodeList[(i+0)*(nx+1) + j]; // bottom

            Cell* negCell = NULL;
            Cell* posCell = NULL;

            if( j != 0  ) negCell = this->CellList[i*nx + (j  - 1)];
            else          negCell = this->CellList[i*nx + (nx - 1)];

            if( j != nx ) posCell = this->CellList[i*nx + (j  - 0)];
            else          posCell = this->CellList[i*nx + (0  - 0)];
            
            bool posAdd = (j != nx);
            bool negAdd = (j != 0) ;

			// create a new interface with the adjacent cells
			Interface* tmpInterface = Interface::createInterface(type,negCell, posCell, negAdd, posAdd, tmpNodes, this->fluidParam, NULL, lengthX);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}

	//=========================================================================
	//=========================================================================
	//						G interface generation
	//=========================================================================
	//=========================================================================
	for (int i = 0; i < ny + 1; i++)        // Y-Direction
	{
		for (int j = 0; j < nx; j++)   // X-Direction
		{
            float2* tmpNodes[2];
            tmpNodes[0] = this->NodeList[i*(nx+1) + (j+0)]; // left
            tmpNodes[1] = this->NodeList[i*(nx+1) + (j+1)]; // right

            Cell* negCell = NULL;
            Cell* posCell = NULL;

            if(i != 0 ) negCell = this->CellList[(i-1)*nx + j];
            if(i != ny) posCell = this->CellList[(i-0)*nx + j];

            BoundaryCondition* currentBC = NULL;
            if( i == 0  ) currentBC = this->BoundaryConditionList[0];
            if( i == ny ) currentBC = this->BoundaryConditionList[1];

			// create a new interface with the adjacent cells
			Interface* tmpInterface = Interface::createInterface(type,negCell, posCell, true, true, tmpNodes, this->fluidParam, currentBC, lengthY);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}

	//=========================================================================
	//=========================================================================
	//						create Ghostcells
	//=========================================================================
	//=========================================================================
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        if( (*i)->isBoundaryInterface() )
        {
            float2* tmpNodes[4];
            float2  normal = (*i)->getScaledNormal();
            
            tmpNodes[1] = (*i)->getNode(0);
            tmpNodes[2] = (*i)->getNode(1);

            tmpNodes[0] = new float2( tmpNodes[1]->x + 2.0*normal.x, tmpNodes[1]->y + 2.0*normal.y );
            tmpNodes[3] = new float2( tmpNodes[2]->x + 2.0*normal.x, tmpNodes[2]->y + 2.0*normal.y );

            this->NodeList.push_back(tmpNodes[0]);
            this->NodeList.push_back(tmpNodes[3]);

            Cell* tmpCell = new Cell(type, tmpNodes, (*i)->getBoundaryCondition(), this->fluidParam);

            tmpCell->addInterface(*i);
            (*i)->addCell(tmpCell);

            this->CellList.push_back(tmpCell);
        }
    }

	//=========================================================================
	//=========================================================================
	//						Compute MinDx
	//=========================================================================
	//=========================================================================
    for(vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        (*i)->computeMinDx();
    }
    
	//=========================================================================
    delete [] CellSpacingsX;
    delete [] NodesX;

    delete [] CellSpacingsY;
    delete [] NodesY;
    
    return;
}

void GKSMesh::generateMiniPatchMesh()
{
    this->lengthX = 1.0;
    this->lengthY = 1.0;

    int nx = 2;
    int ny = 2;

    double eps = 0.1;

    this->NodeList.push_back( new float2( 0.0, 0.0 ) );
    this->NodeList.push_back( new float2( 0.5, 0.0 ) );
    this->NodeList.push_back( new float2( 1.0, 0.0 ) );

    this->NodeList.push_back( new float2( 0.0, 0.5 - eps ) );
    this->NodeList.push_back( new float2( 0.5, 0.5 + eps ) );
    this->NodeList.push_back( new float2( 1.0, 0.5 - eps ) );

    this->NodeList.push_back( new float2( 0.0, 1.0 ) );
    this->NodeList.push_back( new float2( 0.5, 1.0 ) );
    this->NodeList.push_back( new float2( 1.0, 1.0 ) );



	for (int i = 0; i < ny; i++)       // Y-Direction
	{
		for (int j = 0; j < nx; j++)   // X-Direction
		{
            float2* tmpNodes[4];
            tmpNodes[0] = this->NodeList[(i+0)*(nx+1) + j + 0]; // bottom left
            tmpNodes[1] = this->NodeList[(i+0)*(nx+1) + j + 1]; // bottom right
            tmpNodes[2] = this->NodeList[(i+1)*(nx+1) + j + 1]; // top right
            tmpNodes[3] = this->NodeList[(i+1)*(nx+1) + j + 0]; // top left

			Cell* tmpCell = new Cell(compressible, tmpNodes, NULL, this->fluidParam);
			// add interface to list
			this->CellList.push_back(tmpCell);
		}
	}

	for (int i = 0; i < ny; i++)       // Y-Direction
	{
		for (int j = 0; j < nx + 1; j++)    // X-Direction
		{
            float2* tmpNodes[2];
            tmpNodes[0] = this->NodeList[(i+1)*(nx+1) + j]; // top
            tmpNodes[1] = this->NodeList[(i+0)*(nx+1) + j]; // bottom

            Cell* negCell = NULL;
            Cell* posCell = NULL;

            if(j != 0  ) negCell = this->CellList[i*nx + (j  - 1)];
            else         negCell = this->CellList[i*nx + (nx - 1)];

            if(j != nx ) posCell = this->CellList[i*nx + (j  - 0)];
            else         posCell = this->CellList[i*nx + (0  - 0)];
            
            bool posAdd = (j != nx);
            bool negAdd = (j != 0) ;

			// create a new interface with the adjacent cells
			Interface* tmpInterface = Interface::createInterface(compressible,negCell, posCell, negAdd, posAdd, tmpNodes, this->fluidParam, NULL, lengthX);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}
	//=========================================================================
	//=========================================================================
	//						G interface generation
	//=========================================================================
	//=========================================================================
	for (int i = 0; i < ny + 1; i++)        // Y-Direction
	{
		for (int j = 0; j < nx; j++)   // X-Direction
		{
            float2* tmpNodes[2];
            tmpNodes[0] = this->NodeList[i*(nx+1) + (j+0)]; // left
            tmpNodes[1] = this->NodeList[i*(nx+1) + (j+1)]; // right

            Cell* negCell = NULL;
            Cell* posCell = NULL;

            if(i != 0 ) negCell = this->CellList[(i  - 1)*nx + j];
            else        negCell = this->CellList[(ny - 1)*nx + j];
            if(i != ny) posCell = this->CellList[(i  - 0)*nx + j];
            else        posCell = this->CellList[(0  - 0)*nx + j];
            
            bool posAdd = (i != ny);
            bool negAdd = (i != 0) ;

			// create a new interface with the adjacent cells
			Interface* tmpInterface = Interface::createInterface(compressible,negCell, posCell, negAdd, posAdd, tmpNodes, this->fluidParam, NULL, lengthY);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}

	//=========================================================================
	//=========================================================================
	//						Compute MinDx
	//=========================================================================
	//=========================================================================
    for(vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        (*i)->computeMinDx();
    }
}

void GKSMesh::initMeshConstant(double rho, double u, double v, double T)
{
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		(*i)->setValues(rho, u, v, T);
	}
}

void GKSMesh::initMeshLinearTemperature(double rho, double u, double v, double* T)
{
	// Temprature definition
	//    ------------
	//    |   T[1]   |
	//    |          |
	//    |   T[0]   |
	//    ------------
	double interpolatedT;
	float2 center;
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		center = (*i)->getCenter();

        interpolatedT = T[0] + center.y*(T[1] - T[0]) / this->lengthY;

		(*i)->setValues(rho, u, v, interpolatedT);
	}
}

void GKSMesh::initMeshLinear(double* rho, double* u, double* v, double* lambda)
{
	// Temprature definition
	//    ------------
	//    |   Z[1]   |
	//    |          |
	//    |   Z[0]   |
	//    ------------
	double interpolatedRho, interpolatedU, interpolatedV, interpolatedLambda;
	float2 center;
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		center = (*i)->getCenter();

        interpolatedRho    = rho[0]    + center.y*(rho[1]    - rho[0])    / this->lengthY;
        interpolatedU      = u[0]      + center.y*(u[1]      - u[0])      / this->lengthY;
        interpolatedV      = v[0]      + center.y*(v[1]      - v[0])      / this->lengthY;
        //interpolatedLambda = lambda[0] + center.y*(lambda[1] - lambda[0]) / this->lengthY;
        interpolatedLambda = ( lambda[0] ) / ( 1.0 + center.y/this->lengthY * ( lambda[0]/lambda[1] - 1.0 ) );

		(*i)->setValues(interpolatedRho, interpolatedU, interpolatedV, interpolatedLambda);
	}
}

void GKSMesh::initMeshLinearHorizontal(double * rho, double * u, double * v, double * lambda)
{
	// Temprature definition
	//    --------------
	//    |            |
	//    | Z[0]  Z[1] |
	//    |            |
	//    --------------
	double interpolatedRho, interpolatedU, interpolatedV, interpolatedLambda;
	float2 center;
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		center = (*i)->getCenter();

        interpolatedRho    = rho[0]    + center.x*(rho[1]    - rho[0])    / this->lengthX;
        interpolatedU      = u[0]      + center.x*(u[1]      - u[0])      / this->lengthX;
        interpolatedV      = v[0]      + center.x*(v[1]      - v[0])      / this->lengthX;
        //interpolatedLambda = lambda[0] + center.x*(lambda[1] - lambda[0]) / this->lengthX;
        interpolatedLambda = ( lambda[0] ) / ( 1.0 + center.x/this->lengthX * ( lambda[0]/lambda[1] - 1.0 ) );

		(*i)->setValues(interpolatedRho, interpolatedU, interpolatedV, interpolatedLambda);
	}
}

void GKSMesh::initMeshLinearDensity(double * rho, double u, double v, double T)
{
    // Densitsy definition
    //    ------------
    //    |  rho[1]  |
    //    |          |
    //    | rhoT[0]  |
    //    ------------
    double interpolatedRho;
    float2 center;
    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        center = (*i)->getCenter();

        interpolatedRho = rho[0] + center.y*(rho[1] - rho[0]) / this->lengthY;
        //double interpolatedLambda = rho[0] + center.y*( rho[1] - rho[0] ) / this->lengthY;

        //( *i )->setValues(interpolatedRho, u, v, interpolatedLambda);
        (*i)->setValues(interpolatedRho, u, v, T);
    }
}

void GKSMesh::initMeshParabularVelocity(double rho, double u, double v, double T)
{    
    double uValue;
    float2 center;
    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        center = (*i)->getCenter();

        uValue = 4.0 * u * ( center.y - center.y*center.y );

        (*i)->setValues(rho, uValue, v, T);
    }
}

void GKSMesh::initMeshSineVelocity(double rho, double u, double v, double T)
{    
    double uValue;
    float2 center;
    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        center = (*i)->getCenter();

        uValue = u * sin( 2.0 * center.y * M_PI );

        (*i)->setValues(rho, uValue, v, T);
    }
}

void GKSMesh::initMeshAtmospheric(double rho, double u, double v, double lambda, double g)
{
	// Temprature definition
	//    ------------
	//    |   Z[1]   |
	//    |          |
	//    |   Z[0]   |
	//    ------------
	double interpolatedRho;
	float2 center;
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		center = (*i)->getCenter();

        interpolatedRho = rho * exp( -2.0 * g * lambda * center.y );

		(*i)->setValues(interpolatedRho, u, v, lambda);
	}
}

void GKSMesh::addBoundaryCondition( BoundaryConditionType type, 
                                    double rho, double U, double V, double T)
{
    BoundaryCondition* tmp = new BoundaryCondition( type, rho, U, V, T);
    BoundaryConditionList.push_back(tmp);
}

void GKSMesh::addInterfaceBoundaryCondition(double wallVelocity)
{
    InterfaceBC* tmp = new InterfaceBC(wallVelocity);
    this->InterfaceBoundaryConditionsList.push_back(tmp);
}

void GKSMesh::applyBoundaryCondition()
{
    #pragma omp parallel for
    for ( int i = 0; i < CellList.size(); i++ )
    {
        if ( CellList[i]->isGhostCell() )
            CellList[i]->applyBoundaryCondition();
    }
}

void GKSMesh::applyForcing()
{
    #pragma omp parallel for
    for ( int i = 0; i < CellList.size(); i++ )
    {
        if ( ! CellList[i]->isGhostCell() )
            CellList[i]->applyForcing(this->dt);
    }
}

void GKSMesh::computeGlobalTimestep()
{
    this->dt = 1.0e99;
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        if (!((*i)->isGhostCell()))
        {
            this->dt = min( (*i)->getLocalTimestep(), this->dt );
        }
    }
    this->dt *= this->param.CFL;
}

ConservedVariable GKSMesh::getMaxGlobalResidual()
{
    ConservedVariable residual;
    residual.rho  = 0.0;
    residual.rhoU = 0.0;
    residual.rhoV = 0.0;
    residual.rhoE = 0.0;

    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if ( !( ( *i )->isGhostCell() ) )
        {
            residual.rho  = max( residual.rho , (*i)->getLocalResidual().rho  );
            residual.rhoU = max( residual.rhoU, (*i)->getLocalResidual().rhoU );
            residual.rhoV = max( residual.rhoV, (*i)->getLocalResidual().rhoV );
            residual.rhoE = max( residual.rhoE, (*i)->getLocalResidual().rhoE );
        }
    }

    return residual;
}

ConservedVariable GKSMesh::getL2GlobalResidual()
{
    ConservedVariable residual;
    ConservedVariable residualSquare;
    residualSquare.rho  = 0.0;
    residualSquare.rhoU = 0.0;
    residualSquare.rhoV = 0.0;
    residualSquare.rhoE = 0.0;

    ConservedVariable cons;
    cons.rho  = 0.0;
    cons.rhoU = 0.0;
    cons.rhoV = 0.0;
    cons.rhoE = 0.0;

    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if ( !( ( *i )->isGhostCell() ) )
        {
            residualSquare.rho  +=  ( *i )->getLocalResidual().rho  * ( *i )->getLocalResidual().rho;
            residualSquare.rhoU +=  ( *i )->getLocalResidual().rhoU * ( *i )->getLocalResidual().rhoU;
            residualSquare.rhoV +=  ( *i )->getLocalResidual().rhoV * ( *i )->getLocalResidual().rhoV;
            residualSquare.rhoE +=  ( *i )->getLocalResidual().rhoE * ( *i )->getLocalResidual().rhoE;

            cons.rho  +=  ( *i )->getCons().rho  * ( *i )->getCons().rho;
            cons.rhoU +=  ( *i )->getCons().rhoU * ( *i )->getCons().rhoU;
            cons.rhoV +=  ( *i )->getCons().rhoV * ( *i )->getCons().rhoV;
            cons.rhoE +=  ( *i )->getCons().rhoE * ( *i )->getCons().rhoE;
        }
    }

    residual.rho  = sqrt( residualSquare.rho  ) / sqrt( cons.rho  );
    residual.rhoU = sqrt( residualSquare.rhoU ) / sqrt( cons.rhoU );
    residual.rhoV = sqrt( residualSquare.rhoV ) / sqrt( cons.rhoV );
    residual.rhoE = sqrt( residualSquare.rhoE ) / sqrt( cons.rhoE );

    return residual;
}

void GKSMesh::timeStep()
{
    this->iter++;

    if (this->param.verbose) cout << "Iterration: " << this->iter << endl;

    if(this->param.verbose) cout << "  Compute Timestep ..." << endl;
    this->computeGlobalTimestep();
    if (this->param.verbose) cout << "    dt = " << this->dt << endl;

    // ========================================================================

    this->applyForcing();

    this->applyBoundaryCondition();

    // ========================================================================

    if (this->param.verbose) cout << "  Compute Fluxes ..." << endl;

    #pragma omp parallel for
    for ( int i = 0; i < InterfaceList.size(); i++ )
    {
        if ( !InterfaceList[i]->isGhostInterface() )
            InterfaceList[i]->computeFlux(this->dt);
    }

    if (this->param.verbose) cout << "  Update Cells ..." << endl;

    #pragma omp parallel for
    for ( int i = 0; i < CellList.size(); i++ )
    {
        if ( !CellList[i]->isGhostCell() )
            CellList[i]->update(this->dt);
    }

    // ========================================================================

    //this->applyForcing();

    // ========================================================================

    if ( this->param.verbose ) cout << "  Apply Boundary Conditions ..." << endl;
    this->applyBoundaryCondition();

}

void GKSMesh::iterate()
{
    this->time = 0.0;
    this->applyBoundaryCondition();

    ostringstream filename;
    filename << "out/result_0.vtk";
    writeVTKFile(filename.str(), true, false);

    if ( param.fluxOutput == true )
    {
        ostringstream filenameFlux;
        filenameFlux << "out/resultFlux_0.vtk";
        writeVTKFileFlux(filenameFlux.str(), true, false);
    }


    chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

    // ========================================================================
    // ========================================================================
    // ========================================================================
    while (this->iter < this->param.numberOfIterations)
    {
        this->timeStep();
        this->time += this->dt;

        if ( this->iter % this->param.outputInterval == 0 )
        {
            this->dtList.push_back(this->dt);

            ConservedVariable residual = this->getL2GlobalResidual();

            cout << "t = " << this->time << "  \t|  timestep: " << this->iter << endl;
            cout << "r_rho = "  << residual.rho  << "\t ";
            cout << "r_rhoU = " << residual.rhoU << "\t ";
            cout << "r_rhoV = " << residual.rhoV << "\t ";
            cout << "r_rhoE = " << residual.rhoE << "\t ";
            cout << endl;

            this->convergenceHistory.push_back(residual);

            if ( this->isConverged(residual) )
            {
                cout << endl << " ========== Simulation converged! ==========" << endl;
                cout << "Remaining residual change less than " << this->param.convergenceCriterium << endl;
                cout << "Timesteps: " << this->iter << endl;
                cout << "Time: " << this->time << endl;

                ostringstream filename;
                filename << "out/result_" << this->iter << ".vtk";
                writeVTKFile(filename.str(), true, false);

                if ( param.fluxOutput == true )
                {
                    ostringstream filenameFlux;
                    filenameFlux << "out/resultFlux_" << this->iter << ".vtk";
                    writeVTKFileFlux(filenameFlux.str(), true, false);
                }
                break;
            }

        }
        // ========================================================================
        if (this->iter % this->param.outputIntervalVTK == 0)
        {
            ostringstream filename;
            filename << "out/result_" << this->iter << ".vtk";
            writeVTKFile(filename.str(), true, false);

            if ( param.fluxOutput == true )
            {
                ostringstream filenameFlux;
                filenameFlux << "out/resultFlux_" << this->iter << ".vtk";
                writeVTKFileFlux(filenameFlux.str(), true, false);
            }
        }
        // ========================================================================
    }
    // ========================================================================
    // ========================================================================
    // ========================================================================
    
    chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();
    this->computationTime = chrono::duration_cast<chrono::seconds>( endTime - startTime ).count();

    cout << "Time to Solution: " << this->computationTime << " s" << endl;
}

double GKSMesh::getMaxVelocity()
{
    double maxVelocity = 0.0;
    double localVelocity;
    Cell* localCell;
    for ( int i = 0; i < this->CellList.size(); i++ )
    {
        if ( !( this->CellList[i]->isGhostCell() ) )
        {
            localCell = this->CellList[i];
            localVelocity = sqrt( localCell->getPrim().U * localCell->getPrim().U 
                                + localCell->getPrim().V * localCell->getPrim().V );
            maxVelocity = max(maxVelocity, localVelocity);
        }
    }

    return maxVelocity;
}

double GKSMesh::getMaxRe()
{
    return  this->getMaxVelocity()*this->param.L / this->fluidParam.nu;
}

double GKSMesh::getMaxMa()
{
    double maxMa = 0.0;    
    Cell* localCell;
    double localVelocity;
    double localSpeedOfSound;

    for ( int i = 0; i < this->CellList.size(); i++ )
    {
        if ( !( this->CellList[i]->isGhostCell() ) )
        {
            localCell = this->CellList[i];
            localVelocity = sqrt( localCell->getPrim().U * localCell->getPrim().U
                                + localCell->getPrim().V * localCell->getPrim().V );
            localSpeedOfSound = sqrt(1.0 / ( 2.0 * localCell->getPrim().L ));
            maxMa = max(maxMa, localVelocity/localSpeedOfSound);
        }
    }

    return maxMa;
}

bool GKSMesh::isConverged(ConservedVariable residual)
{
    bool flag = true;

    flag = flag && ( residual.rho  < this->param.convergenceCriterium[0] );
    flag = flag && ( residual.rhoU < this->param.convergenceCriterium[1] );
    flag = flag && ( residual.rhoV < this->param.convergenceCriterium[2] );
    flag = flag && ( residual.rhoE < this->param.convergenceCriterium[3] );

    return flag;
}

string GKSMesh::toString()
{
	ostringstream tmp;
	tmp << "The Mesh has following interfaces:\n";

	for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
	{
        //if( ! (*i)->isGhostInterface() )
		tmp << (*i)->toString();
	}

	return tmp.str();
}

string GKSMesh::cellValuesToString()
{
    ostringstream tmp;
    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        if( ! (*i)->isGhostCell()  )
            tmp << (*i)->valuesToString() << "\n";
    }
    return tmp.str();
}

void GKSMesh::writeOverviewFile(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.precision(15);
    file.open(filename.c_str());

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    file << " ========== Fluid Parameters ==========";
    file << endl;
    file << "nu =\t " << this->fluidParam.nu << endl;
    file << "K  =\t " << this->fluidParam.K << endl;
    file << "R  =\t " << this->fluidParam.R << endl;
    file << "Fx =\t " << this->fluidParam.Force.x << endl;
    file << "Fy =\t " << this->fluidParam.Force.y << endl;
    file << endl;

    file << " ========== Simulation Parameters ==========";
    file << endl;
    file << "Max Number of Iteratios:             " << this->param.numberOfIterations << endl;
    file << "VTK-File Output Intervat:            " << this->param.outputIntervalVTK << endl;
    file << "Convergence History Output Interval: " << this->param.outputInterval << endl;
    file << "Convergence Criterium:               ( " << this->param.convergenceCriterium[0] << ", "
                                                      << this->param.convergenceCriterium[1] << ", "
                                                      << this->param.convergenceCriterium[2] << ", "
                                                      << this->param.convergenceCriterium[3] << " )" << endl;
    file << "CFL =\t" << this->param.CFL << endl;
    file << endl;

    file << " ========== Flow Characteristics ==========";
    file << endl;
    file << "Umax = " << this->getMaxVelocity() << " m/s" << endl;
    file << "Re   = " << this->getMaxRe() << endl;
    file << "Ma   = " << this->getMaxMa() << endl;
    file << endl;

    file << " ========== Simulation Results ==========";
    file << endl;
    if ( this->isConverged(*(convergenceHistory.end()-1)) )
    {
        file << "Simulation Converged!" << endl;
        file << "Time steps: " << this->iter << endl;
    }
    else
    {
        file << "Simulation did not converge" << endl;
    }
    file << "Final Residuals:" << endl;
    file << "r_rho  = " << ( *(convergenceHistory.end()-1) ).rho << "\t ";
    file << "r_rhoU = " << ( *(convergenceHistory.end()-1) ).rhoU << "\t ";
    file << "r_rhoV = " << ( *(convergenceHistory.end()-1) ).rhoV << "\t ";
    file << "r_rhoE = " << ( *(convergenceHistory.end()-1) ).rhoE << "\t ";
    file << endl;
    file << endl;
    file << "Real Time Simulated : " << this->time << " s" << endl;
    file << "Time to solution:     " << this->computationTime << " s" << endl;
    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeVTKFile(string filename, bool data, bool BC)
{
    cout << "Wrinting file " << filename << " ... ";
	// open file stream
	ofstream file;
    file.precision(15);
	file.open(filename.c_str());

	if (!file.is_open()) {
		cout << " File cound not be opened.\n\nERROR!\n\n\n";
		return;
	}

    this->writeCellGeometry(file);

    if (data) this->writeCellData(file);

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeVTKFileFlux(string filename, bool data, bool BC)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.precision(15);
    file.open(filename.c_str());

    if (!file.is_open()) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    this->writeInterfaceGeometry(file);

    if (data) this->writeInterfaceData(file);

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeTimeSteps(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.precision(15);
    file.open(filename.c_str());

    if (!file.is_open()) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for (vector<double>::iterator i = this->dtList.begin(); i != this->dtList.end(); ++i)
        file << (*i) << "\n";

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeVelocityProfile(string filename, double x)
{

    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.precision(15);
    file.open(filename.c_str());

    if (!file.is_open()) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        if ( !( *i )->isGhostCell() )
        {
            double xMax = -1e99;
            double xMin =  1e99;
            for(int j = 0; j < 4; j++)
            {
                if( ( *i )->getNode(j).x > xMax) xMax = ( *i )->getNode(j).x;
                if( ( *i )->getNode(j).x < xMin) xMin = ( *i )->getNode(j).x;
            }

            // check wether the profile location x is located in this cell
            if ( x <= xMax && x >= xMin )
            {
                file << ( *i )->getCenter().y << " " << ( *i )->getPrim().U << " " << ( *i )->getPrim().rho << "\n";
            }
        }
    }

    file.close();

    cout << "done!" << endl;

}

void GKSMesh::writeResultFields(string filename)
{

    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.precision(15);
    file.open(filename.c_str());

    if (!file.is_open()) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        file << ( *i )->getCenter().x << " ";
        file << ( *i )->getCenter().y << " ";
        file << ( *i )->getPrim().rho << " ";
        file << ( *i )->getPrim().U << " ";
        file << ( *i )->getPrim().V << " ";
        file << 1.0 / ( 2.0 * ( *i )->getPrim().L * this->fluidParam.R ) << " ";
        file << ( *i )->isGhostCell() << endl;
    }

    file.close();

    cout << "done!" << endl;

}

void GKSMesh::writeTemperatureProfile(string filename, double x)
{

    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.precision(15);
    file.open(filename.c_str());

    if (!file.is_open()) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        if ( !( *i )->isGhostCell() )
        {
            double xMax = -1e99;
            double xMin =  1e99;
            for(int j = 0; j < 4; j++)
            {
                if( ( *i )->getNode(j).x > xMax) xMax = ( *i )->getNode(j).x;
                if( ( *i )->getNode(j).x < xMin) xMin = ( *i )->getNode(j).x;
            }

            // check wether the profile location x is located in this cell
            if ( x <= xMax && x >= xMin )
            {
                file << ( *i )->getCenter().y << " " << 1.0 / ( 2.0 * this->fluidParam.R * ( *i )->getPrim().L ) << "\n";
            }
        }
    }

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeMeshAsText(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.open(filename.c_str());

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    file << " ========== Cells: ========== \n\n";

    for ( vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i )
    {
        file << (*i)->toString();
        if ( ( *i )->isGhostCell() )
            file << ", Ghostcell";
        file << "\n";
    }

    file << "\n\n ========== Interfaces: ========== \n\n";

    for ( vector<Interface*>::iterator i = this->InterfaceList.begin(); i != this->InterfaceList.end(); ++i )
    {
        file << ( *i )->toString();
    }

    file.close();
}

void GKSMesh::writeVelocityU(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.open(filename.c_str());

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        //if ( !( *i )->isGhostCell() )
            file << ( *i )->getPrim().U << endl;
    }

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeVelocityV(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.open(filename.c_str());

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        //if( !(*i)->isGhostCell() )
            file << ( *i )->getPrim().V << endl;
    }

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeTemperature(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.open(filename.c_str());

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        //if( !(*i)->isGhostCell() )
            file << 1.0 / (2.0 * this->fluidParam.R * ( *i )->getPrim().L ) << endl;
    }

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeDensity(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.open(filename.c_str());

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        //if( !(*i)->isGhostCell() )
            file << ( *i )->getPrim().rho << endl;
    }

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeConvergenceHistory(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.precision(15);
    file.open(filename.c_str());

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }

    for ( int i = 0; i < this->convergenceHistory.size(); i++ )
    {
        file << convergenceHistory[i].rho  << "\t";
        file << convergenceHistory[i].rhoU << "\t";
        file << convergenceHistory[i].rhoV << "\t";
        file << convergenceHistory[i].rhoE << "\t";
        file << endl;
    }

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeCellGeometry(ofstream& file)
{

    // write VTK Header
    file << "# vtk DataFile Version 1.0\n";
    file << "by Stephan Lenz\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // write nodes
    //( one dummy node with the ID 0 must be written )
    file << "POINTS " << 4 * this->CellList.size() + 1 << " double\n";
    file << "0.0 0.0 0.0 \n";

    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->writeNodes();
    }

    // write elements
    file << "CELLS " << this->CellList.size() << " " << 5 * this->CellList.size() << endl;
    for (int i = 0; i < this->CellList.size(); ++i)
    {
        file << 4 << " " << i * 4 + 1
            << " " << i * 4 + 2
            << " " << i * 4 + 3
            << " " << i * 4 + 4 << endl;
    }

    // write element tyes( 9 = quad element )
    file << "CELL_TYPES " << this->CellList.size() << endl;
    for (int i = 0; i < this->CellList.size(); ++i)
    {
        file << "9" << endl;
    }
}

void GKSMesh::writeInterfaceGeometry(ofstream& file)
{

    // write VTK Header
    file << "# vtk DataFile Version 1.0\n";
    file << "by Stephan Lenz\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // write nodes
    file << "POINTS " << this->InterfaceList.size()<< " double\n";

    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << (*i)->writeCenter();
    }

    // write elements
    file << "CELLS " << this->InterfaceList.size() << " " << 2 * this->InterfaceList.size() << endl;
    for (int i = 0; i < this->InterfaceList.size(); ++i)
    {
        file << 1 << " " << i << endl;
    }

    // write element tyes( 9 = quad element )
    file << "CELL_TYPES " << this->InterfaceList.size() << endl;
    for (int i = 0; i < this->InterfaceList.size(); ++i)
    {
        file << "1" << endl;
    }
}

void GKSMesh::writeCellData(ofstream& file)
{
    int numberOfFields = 10;
    if ( this->param.resOutput )
        numberOfFields += 4;
    
    // write cell data ( ID and stress )
    file << "CELL_DATA " << this->CellList.size() << endl;
    file << "FIELD Lable " << numberOfFields << "\n";

    // ================================================================================================================

    file << "rho 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getPrim().rho << endl;
    }

    file << "U 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getPrim().U << endl;
    }

    file << "V 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getPrim().V << endl;
    }

    file << "Lambda 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getPrim().L << endl;
    }

    // ================================================================================================================

    file << "GhostCell 1 " << this->CellList.size() << " int\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if ( ( *i )->isGhostCell() )
            file << 1 << endl;
        else
            file << 0 << endl;
    }

    // ================================================================================================================

    file << "rhoU 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getCons().rhoU << endl;
    }

    file << "rhoV 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getCons().rhoV << endl;
    }

    file << "rhoE 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getCons().rhoE << endl;
    }

    // ================================================================================================================

    file << "p 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getPrim().rho / ( 2.0 * ( *i )->getPrim().L ) << endl;
    }

    file << "T 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << 1.0 / ( 2.0 * this->fluidParam.R * ( *i )->getPrim().L ) << endl;
    }

    // ================================================================================================================

    if ( this->param.resOutput )
    {
        file << "res_rho 1 " << this->CellList.size() << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            file << ( *i )->getLocalResidual().rho << endl;
        }

        file << "res_rhoU 1 " << this->CellList.size() << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            file << ( *i )->getLocalResidual().rhoU << endl;
        }

        file << "res_rhoV 1 " << this->CellList.size() << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            file << ( *i )->getLocalResidual().rhoV << endl;
        }

        file << "res_rhoE 1 " << this->CellList.size() << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            file << ( *i )->getLocalResidual().rhoE << endl;
        }
    }
    // ================================================================================================================
    for(int j = 0; j < 4; j++)
    {
        file << "VECTORS Connectivity" << j << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            file << (*i)->getConnectivity(j).x << " " << (*i)->getConnectivity(j).y << " 0.0" << endl;
        }
    }
}

void GKSMesh::writeInterfaceData(ofstream & file)
{
    // write cell data ( ID and stress )
    file << "POINT_DATA " << this->InterfaceList.size() << endl;
    file << "FIELD Lable 8\n";

    file << "F_rho 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        //file << (*i)->getTimeIntegratedFlux().rho << endl;
        file << ( *i )->getFluxDensity().rho << endl;
    }

    file << "F_rhoU 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        //file << (*i)->getTimeIntegratedFlux().rhoU << endl;
        file << ( *i )->getFluxDensity().rhoU << endl;
    }

    file << "F_rhoV 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        //file << (*i)->getTimeIntegratedFlux().rhoV << endl;
        file << ( *i )->getFluxDensity().rhoV << endl;
    }

    file << "F_rhoE 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        //file << (*i)->getTimeIntegratedFlux().rhoE << endl;
        file << ( *i )->getFluxDensity().rhoE << endl;
    }

    file << "F_rho_t 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << (*i)->getTimeIntegratedFlux().rho << endl;
    }

    file << "F_rhoU_t 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << (*i)->getTimeIntegratedFlux().rhoU << endl;
    }

    file << "F_rhoV_t 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << (*i)->getTimeIntegratedFlux().rhoV << endl;
    }

    file << "F_rhoE_t 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << (*i)->getTimeIntegratedFlux().rhoE << endl;
    }

    //file << "GhostInterface 1 " << this->InterfaceList.size() << " int\n";
    //for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    //{
    //    if ((*i)->isGhostInterface())
    //        file << 1 << endl;
    //    else
    //        file << 0 << endl;
    //}

    file << "VECTORS normal double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << (*i)->getNormal().x << " " << (*i)->getNormal().y << " 0.0" << endl;
    }

    file << "VECTORS posCell double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << (*i)->getPosConnectivity().x << " " << (*i)->getPosConnectivity().y << " 0.0" << endl;
    }

    file << "VECTORS negCell double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << (*i)->getNegConnectivity().x << " " << (*i)->getNegConnectivity().y << " 0.0" << endl;
    }
}

