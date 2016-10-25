// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      GKSMesh.cpp
//
//      Function:
//          Generation Storage of mesh
//          Control of simulation
//          Data Analysis
//          File Output
//
// ============================================================================

#include "GKSMesh.h"
#include "Cell.h"
#include "CompressibleInterface.h"
#include "Types.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>    // min()
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

// ============================================================================
//              deletes all dynamically allocated memory
// ============================================================================
GKSMesh::~GKSMesh()
{
    // ========================================================================
    //              delete Node Memory
    // ========================================================================
    for( int i = 0; i < this->NodeList.size(); i++ )
        delete this->NodeList[i];
    
    // ========================================================================
    //              delete Cell Memory
    // ========================================================================
    for( int i = 0; i < this->CellList.size(); i++ )
        delete this->CellList[i];
    
    // ========================================================================
    //              delete Node Memory
    // ========================================================================
    for( int i = 0; i < this->InterfaceList.size(); i++ )
        delete this->InterfaceList[i];
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          Mesh generation
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method generates a structured mesh on a rectangular domain.
//      
//      Parameters:
//          type:       either compressible or incompressible
//          lengthX:    length of the domain in x-direction
//          lengthY:    length of the domain in y-direction
//          nx:         number of cells in x-direction
//          ny:         number of cells in y-direction
//          gradingX:   ratio of the inner cells to the outer cells in x direction
//          gradingY:   ratio of the inner cells to the outer cells in y direction
//          
// ============================================================================
void GKSMesh::generateRectMeshGraded(InterfaceType type, double lengthX, double lengthY, int nx, int ny, double gradingX, double gradingY)
{

    this->lengthX = lengthX;
    this->lengthY = lengthY;

    // make nx and ny an even number for two sided grading
    nx = 2 * (nx/2);
    ny = 2 * (ny/2);
    
    // compute expansion coefficients eta for cell size (dx_n+1 = eta * dx_n)
    double etaX = pow( gradingX, 1.0 / (nx/2 - 1) );
    double etaY = pow( gradingY, 1.0 / (ny/2 - 1) );

    // ========================================================================
    //      compute size of the inner cells by solution of partial geometric series
    // ========================================================================
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
    // ========================================================================

    //=========================================================================
    //		Computate the cell sizes
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
    
    //=========================================================================
    //		Computate the node coordinates
    //=========================================================================
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
	//=========================================================================
	//		Node generation
	//=========================================================================
	//=========================================================================
    int currentNodeID = 1;

    for (int i = 0; i < ny + 1; i++)       // Y-Direction
	{
		for (int j = 0; j < nx + 1; j++)   // X-Direction
		{
            // ===== No Distortion =====================
            float2* tmpNode = new float2( NodesX[j], NodesY[i] );
            
            tmpNode->ID = currentNodeID++;

            // ================================================================
            //          in this section several mesh distortions can be
            //          applied to the Mesh
            // ================================================================

            // ===== Parallelogram =====================
            //tmpNode->y += NodesX[j] / this->lengthX * 0.5;

            // ===== internal x-Distortion =============
            //tmpNode->x += - 2.0 * (NodesY[i] - this->lengthY)*NodesY[i]*(NodesX[j] - this->lengthX)*NodesX[j];

            // ===== x-Distortion (parallel) ===========
            //tmpNode->x += - 0.4 * this->lengthX * sin( NodesX[j] / this->lengthX * M_PI );

            // ===== internal y-Distortion =============
            //tmpNode->y -= - 2.0 * (NodesX[j] - this->lengthX)*NodesX[j]*(NodesY[i] - this->lengthY)*NodesY[i];
            //tmpNode->y += 0.2 * sin( NodesX[j] * M_PI ) * sin( NodesY[i] * M_PI );

            // ===== internal y-Distortion (symmetric) =
            //tmpNode->y += 0.05 * cos( NodesX[j] * 2.0 * M_PI / this->lengthX ) * sin( NodesY[i] * M_PI/this->lengthY );

            // ===== y-Distortion (parallel) ===========
            //tmpNode->y += - 0.4 * this->lengthY * sin( NodesY[i] / this->lengthY * M_PI );

            // ===== internal parabular Distortion =====
            //tmpNode->y += - 0.25 * (NodesX[j] - this->lengthX)*NodesX[j] * sin( (NodesY[i] - 0.5*this->lengthY) * 2.0 * M_PI/this->lengthY );

            // ===== internal sine Distortion (x) ======
            //tmpNode->x -=  0.05 * sin( NodesY[i] * 2.0 * M_PI/this->lengthY ) * sin( (NodesX[j] - 0.5*this->lengthX) * 2.0 * M_PI/this->lengthX );

            // ===== internal sine Distortion (y) ======
            //tmpNode->y -= 0.05 * sin( NodesX[j] * 2.0 * M_PI/this->lengthX ) * sin( ( NodesY[i] - 0.5*this->lengthY ) * 2.0 * M_PI/this->lengthY );

            // ===== internal linear Distortion ========
            //if( i != 0 && i != ny && j != 0 && j != nx ) tmpNode->y += 0.1;

			this->NodeList.push_back(tmpNode);
		}
	}
	//=========================================================================
    
	//=========================================================================
	//=========================================================================
	//		Node Rotation - can be used to rotate the whole mesh
	//=========================================================================
	//=========================================================================
    double angle = 0.0;         // rotation angle in arc

    for (vector<float2*>::iterator i = NodeList.begin(); i != NodeList.end(); ++i)
	{
        double x = (*i)->x;
        double y = (*i)->y;
        (*i)->x = cos(angle) * x - sin(angle) * y;
        (*i)->y = sin(angle) * x + cos(angle) * y;
	}
	//=========================================================================

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

			// add pointer of the cell to CellList
			this->CellList.push_back(tmpCell);
		}
	}
	//=========================================================================

	//=========================================================================
	//=========================================================================
	//                      F interface generation
	//=========================================================================
	//=========================================================================
    for (int i = 0; i < ny; i++)            // Y-Direction
    {
        for (int j = 0; j < nx + 1; j++)    // X-Direction
        {
            float2* tmpNodes[2];
            tmpNodes[0] = this->NodeList[(i+1)*(nx+1) + j]; // top
            tmpNodes[1] = this->NodeList[(i+0)*(nx+1) + j]; // bottom

            Cell* negCell = NULL;
            Cell* posCell = NULL;

            // ================================================================
            // find adjacent cells for this interface
            //      in case of periodic boundaries this is the cell on the 
            //      opposite side of the domain
            // ================================================================
            if( j != 0  ) 
                negCell = this->CellList[i*nx + (j  - 1)];
            else if( this->BoundaryConditionList[0]->getType() == periodic || this->BoundaryConditionList[0]->getType() == periodicGhost )
                negCell = this->CellList[i*nx + (nx - 1)];

            if( j != nx ) 
                posCell = this->CellList[i*nx + (j  - 0)];
            else if( this->BoundaryConditionList[2]->getType() == periodic || this->BoundaryConditionList[2]->getType() == periodicGhost )       
                posCell = this->CellList[i*nx + (0  - 0)];
            
            // ================================================================
            //      in case of periodic boundaries without ghost cells, the
            //      interface does not need to know wether it is on the
            //      boundary or not
            // ================================================================
            BoundaryCondition* currentBC = NULL;
            if( j == 0  && this->BoundaryConditionList[0]->getType() != periodic )
                currentBC = this->BoundaryConditionList[0];

            if( j == nx && this->BoundaryConditionList[2]->getType() != periodic )
                currentBC = this->BoundaryConditionList[2];
            
            // ================================================================
            //      these flags tell, if this interface shall be added to the 
            //      pos/neg cells or not
            // ================================================================
            bool posAdd = (j != nx);
            bool negAdd = (j != 0) ;

			// create a new interface with the adjacent cells
			Interface* tmpInterface = Interface::createInterface(type,negCell, posCell, negAdd, posAdd, tmpNodes, this->fluidParam, currentBC, lengthX, lengthY);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}
	//=========================================================================

	//=========================================================================
	//=========================================================================
	//                      G interface generation
	//=========================================================================
	//=========================================================================
	for (int i = 0; i < ny + 1; i++)   // Y-Direction
	{
		for (int j = 0; j < nx; j++)   // X-Direction
		{
            float2* tmpNodes[2];
            tmpNodes[0] = this->NodeList[i*(nx+1) + (j+0)]; // left
            tmpNodes[1] = this->NodeList[i*(nx+1) + (j+1)]; // right

            Cell* negCell = NULL;
            Cell* posCell = NULL;

            // ================================================================
            // find adjacent cells for this interface
            //      in case of periodic boundaries this is the cell on the 
            //      opposite side of the domain
            // ================================================================
            if( i != 0  ) 
                negCell = this->CellList[(i -1)*nx + j];
            else if( this->BoundaryConditionList[1]->getType() == periodic || this->BoundaryConditionList[1]->getType() == periodicGhost )
                negCell = this->CellList[(ny-1)*nx + j];

            if( i != ny ) 
                posCell = this->CellList[(i -0)*nx + j];
            else if( this->BoundaryConditionList[3]->getType() == periodic || this->BoundaryConditionList[3]->getType() == periodicGhost )       
                posCell = this->CellList[(0 -0)*nx + j];
            
            // ================================================================
            //      in case of periodic boundaries without ghost cells, the
            //      interface does not need to know wether it is on the
            //      boundary or not
            // ================================================================
            BoundaryCondition* currentBC = NULL;
            if( i == 0  && this->BoundaryConditionList[1]->getType() != periodic ) 
                currentBC = this->BoundaryConditionList[1];

            if( i == ny && this->BoundaryConditionList[3]->getType() != periodic ) 
                currentBC = this->BoundaryConditionList[3];
            
            // ================================================================
            //      these flags tell, if this interface shall be added to the 
            //      pos/neg cells or not
            // ================================================================
            bool posAdd = (i != ny);
            bool negAdd = (i != 0) ;

			// create a new interface with the adjacent cells
			Interface* tmpInterface = Interface::createInterface(type,negCell, posCell, negAdd, posAdd, tmpNodes, this->fluidParam, currentBC, lengthX, lengthY);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}
	//=========================================================================

	//=========================================================================
	//=========================================================================
	//						create Ghostcells
    //                          Ghost cells are created for all Interfaces
    //                          that have a boundary condition pointer
	//=========================================================================
	//=========================================================================
    vector<Interface*> tmpInterfaceList;

    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        if( (*i)->isBoundaryInterface() )
        {
            Cell* periodicCell = NULL;
            float2* tmpNodes[4];
            
            tmpNodes[1] = (*i)->getNode(0);
            tmpNodes[2] = (*i)->getNode(1);

            if( (*i)->getBoundaryCondition()->getType() == periodicGhost)
            {
                periodicCell = (*i)->getPeriodicCell();
            }
            
            // ================================================================
            //      create the new nodes such, that the cell center of the 
            //      ghost cell has the same distance to the interface as the 
            //      cell center of the adjacent cell in the domain
            // ================================================================
            float2  normal = (*i)->getScaledNormal();
            tmpNodes[0] = new float2( tmpNodes[1]->x + 2.0*normal.x, tmpNodes[1]->y + 2.0*normal.y );
            tmpNodes[3] = new float2( tmpNodes[2]->x + 2.0*normal.x, tmpNodes[2]->y + 2.0*normal.y );

            tmpNodes[0]->ID = 9999999; // 9.999.999
            tmpNodes[3]->ID = 9999999; // 9.999.999

            this->NodeList.push_back(tmpNodes[0]);
            this->NodeList.push_back(tmpNodes[3]);

            Cell* tmpCell = new Cell(type, tmpNodes, (*i)->getBoundaryCondition(), this->fluidParam);

            tmpCell->addInterface(*i);
            (*i)->addCell(tmpCell);

            // ================================================================
            //      In the case of periodic boundaries with ghost cells, the 
            //      ghost cell is connected to its counterpart in the domain
            //      on the opposite side of the domain by an interface. This
            //      Interface is created here.
            // ================================================================
            if( (*i)->getBoundaryCondition()->getType() == periodicGhost)
            {
                float2* tmpNodesInterface[2];
                tmpNodesInterface[0] = tmpNodes[0];
                tmpNodesInterface[1] = tmpNodes[3];

                Interface* tmpInterface = Interface::createInterface(type, tmpCell, periodicCell, true, true, tmpNodesInterface, fluidParam, (*i)->getBoundaryCondition(), lengthX, lengthY);

                tmpInterfaceList.push_back(tmpInterface);
            }

            this->CellList.push_back(tmpCell);
        }
    }

    this->InterfaceList.insert( this->InterfaceList.end(), tmpInterfaceList.begin(), tmpInterfaceList.end() );
    tmpInterfaceList.clear();
	//=========================================================================

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

	//=========================================================================
	//=========================================================================
	//						Compute least square coefficients
	//=========================================================================
	//=========================================================================
    for(vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        if( !(*i)->isGhostCell() )
        {
            (*i)->computeLeastSquareCoefficients();
        }
    }
    
	//=========================================================================
    delete [] CellSpacingsX;
    delete [] NodesX;

    delete [] CellSpacingsY;
    delete [] NodesY;

	return;
}

// ============================================================================
//      This method adds a Boundary Condition to the mesh. This must be done
//      before the mesh generation.
//
//      Parameters:
//          type:   wall, isothermalWall, periodic, periodicGhost
//          rho:    density on the boundary         (currently not used)
//          U:      x-velocity on the boundary      (only used for wall and isothermalWall)
//          V:      y-velocity on the boundary      (only used for wall and isothermalWall)
//          T:      Temperature on the boundary     (only used for isothermalWall)
// ============================================================================
void GKSMesh::addBoundaryCondition( BoundaryConditionType type, double rho, double U, double V, double T)
{
    BoundaryCondition* tmp = new BoundaryCondition( type, rho, U, V, T);
    BoundaryConditionList.push_back(tmp);
}

// ============================================================================
//      Interface Boundary Conditions are currently not supported
// ============================================================================
void GKSMesh::addInterfaceBoundaryCondition(double wallVelocity)
{
    InterfaceBC* tmp = new InterfaceBC(wallVelocity);
    this->InterfaceBoundaryConditionsList.push_back(tmp);
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          Mesh initialization
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method initializes the mesh with constant value every where.
//
//      Parameters:
//          rho:    density
//          U:      x-velocity
//          V:      Y-velocity
//          T:      Temperature
// ============================================================================
void GKSMesh::initMeshConstant(double rho, double U, double V, double T)
{
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		(*i)->setValues(rho, U, V, T);
	}
}

// ============================================================================
//      This method initializes the mesh with linear distribution in
//      y-direction.
//          Value definition, Z = (rho,U,V,T)
//             ------------
//             |   Z[1]   |
//             |          |
//             |   Z[0]   |
//             ------------
//
//      Parameters:
//          rho:    density
//          U:      x-velocity
//          V:      Y-velocity
//          T:      Temperature
// ============================================================================
void GKSMesh::initMeshLinear(double* rho, double* U, double* V, double* lambda)
{
	double interpolatedRho, interpolatedU, interpolatedV, interpolatedLambda;
	float2 center;
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		center = (*i)->getCenter();

        interpolatedRho    = rho[0]    + center.y*(rho[1]    - rho[0])    / this->lengthY;
        interpolatedU      = U[0]      + center.y*(U[1]      - U[0])      / this->lengthY;
        interpolatedV      = V[0]      + center.y*(V[1]      - V[0])      / this->lengthY;
        interpolatedLambda = ( lambda[0] ) / ( 1.0 + center.y/this->lengthY * ( lambda[0]/lambda[1] - 1.0 ) );

		(*i)->setValues(interpolatedRho, interpolatedU, interpolatedV, interpolatedLambda);
	}
}

// ============================================================================
//      This method initializes the mesh with linear distribution in
//      x-direction.
//          Value definition, Z = (rho,U,V,T)
//             --------------
//             |            |
//             | Z[0]  Z[1] |
//             |            |
//             --------------
//
//      Parameters:
//          rho:    density
//          U:      x-velocity
//          V:      Y-velocity
//          T:      Temperature
// ============================================================================
void GKSMesh::initMeshLinearHorizontal(double * rho, double * U, double * V, double * lambda)
{  
	double interpolatedRho, interpolatedU, interpolatedV, interpolatedLambda;
	float2 center;
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		center = (*i)->getCenter();

        interpolatedRho    = rho[0]    + center.x*(rho[1]    - rho[0])    / this->lengthX;
        interpolatedU      = U[0]      + center.x*(U[1]      - U[0])      / this->lengthX;
        interpolatedV      = V[0]      + center.x*(V[1]      - V[0])      / this->lengthX;
        interpolatedLambda = ( lambda[0] ) / ( 1.0 + center.x/this->lengthX * ( lambda[0]/lambda[1] - 1.0 ) );

		(*i)->setValues(interpolatedRho, interpolatedU, interpolatedV, interpolatedLambda);
	}
}

// ============================================================================
//      This method initializes the mesh with parabolic x-velocity with U in 
//      the middle and 0 at top and bottom. The other quantities are constant.
//
//      Parameters:
//          rho:    density
//          U:      x-velocity
//          V:      Y-velocity
//          T:      Temperature
// ============================================================================
void GKSMesh::initMeshParabularVelocity(double rho, double U, double V, double T)
{    
    double uValue;
    float2 center;
    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        center = (*i)->getCenter();

        uValue = 4.0 * U * ( center.y - center.y*center.y );

        (*i)->setValues(rho, uValue, V, T);
    }
}


// ============================================================================
//      This method initializes the mesh with sinusodial x-velocity with U 
//      being the maximal x velocity. The other quantities are constant.
//
//      Parameters:
//          rho:    density
//          U:      x-velocity
//          V:      Y-velocity
//          T:      Temperature
// ============================================================================
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

// ============================================================================
//      This method initializes the mesh with Such that it is stable under
//      constant external forcing (gravity).
//          Density definition
//             -----------
//             |         |
//             |         |
//             |   rho   |
//             -----------
//
//      Parameters:
//          rho:    density at the bottom
//          U:      x-velocity
//          V:      Y-velocity
//          T:      Temperature
// ============================================================================
void GKSMesh::initMeshAtmospheric(double rho, double u, double v, double lambda, double g)
{
	double interpolatedRho;
	float2 center;
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		center = (*i)->getCenter();

        interpolatedRho = rho * exp( -2.0 * g * lambda * center.y );

		(*i)->setValues(interpolatedRho, u, v, lambda);
	}
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          Simulation Control
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method perform the iteration and controls the solution process. 
//      It also call the file output and checks the convergence.
// ============================================================================
void GKSMesh::iterate()
{
    this->iter = 0;
    this->time = 0.0;

    this->timeList.push_back(this->time);

    this->writeOutputFiles();

    chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

    // ========================================================================
    // ========================================================================
    //              Time Loop
    // ========================================================================
    // ========================================================================
    while (this->iter < this->param.numberOfIterations)
    {
        // ====================================================================
        this->iter++;
        // ====================================================================

        // ====================================================================
        //          Perform one timestep
        // ====================================================================
        this->timeStep();
        // ====================================================================

        this->time += this->dt;
        
        // ====================================================================
        if ( this->iter % this->param.outputInterval == 0 )
        {
            this->dtList.push_back(this->dt);
            this->timeList.push_back(this->time);

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

                this->writeOutputFiles();

                break;
            }
        }
        // ====================================================================

        // ====================================================================
        if (this->iter % this->param.outputIntervalVTK == 0)
        {
            this->writeOutputFiles();
        }
        // ====================================================================
    }
    // ========================================================================
    // ========================================================================
    //              End of Time Loop
    // ========================================================================
    // ========================================================================
    
    chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();
    this->computationTime = chrono::duration_cast<chrono::seconds>( endTime - startTime ).count();

    cout << "Time to Solution: " << this->computationTime << " s" << endl;
}

// ============================================================================
//      This method performes on time step of the GKS.
//      All computation is taking place with in this method.
// ============================================================================
void GKSMesh::timeStep()
{
    this->computeGlobalTimestep();

    this->applyForcing();

    this->applyBoundaryCondition();

    this->computeLeastSquareGradients();

    this->applyBoundaryCondition();

    this->computeFluxes();

    this->updateCells();

    this->applyBoundaryCondition();
}

// ============================================================================
//      This method computes the global CFL time step by finding the minimum
//      of all local time steps and multiplying it with the CFL number.
// ============================================================================
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

// ============================================================================
//      This method applies the forcing to all cells
// ============================================================================
void GKSMesh::applyForcing()
{
    #pragma omp parallel for
    for ( int i = 0; i < CellList.size(); i++ )
    {
        if ( ! CellList[i]->isGhostCell() )
            CellList[i]->applyForcing(this->dt);
    }
}

// ============================================================================
//      This method sets the values of the ghost cells
// ============================================================================
void GKSMesh::applyBoundaryCondition()
{
    #pragma omp parallel for
    for ( int i = 0; i < CellList.size(); i++ )
    {
        if ( CellList[i]->isGhostCell() )
            CellList[i]->applyBoundaryCondition();
    }
}

// ============================================================================
//      This method computes the gradients of the conserved variables in the 
//      Cells with the least suqare method. These gradients are not used, if
//      Finite differences are applied at the Interfaces.
// ============================================================================
void GKSMesh::computeLeastSquareGradients()
{
    #pragma omp parallel for
    for ( int i = 0; i < CellList.size(); i++ )
    {
        // TODO: Right now the Gradients in the periodic ghost cells are 
        // computed by finite differences and not by least square
        //if ( !CellList[i]->isGhostCell() )
            CellList[i]->computeGradients();
    }
}

// ============================================================================
//      This method computes the fluxes over all Interfaces
// ============================================================================
void GKSMesh::computeFluxes()
{
    #pragma omp parallel for
    for ( int i = 0; i < InterfaceList.size(); i++ )
    {
        InterfaceList[i]->computeFlux(this->dt);
    }
}

// ============================================================================
//      This method updates the cell averaged conserved variables from the fluxes.
// ============================================================================
void GKSMesh::updateCells()
{
    #pragma omp parallel for
    for ( int i = 0; i < CellList.size(); i++ )
    {
        if ( !CellList[i]->isGhostCell() )
            CellList[i]->update(this->dt);
    }
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          Data Analysis
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method returns the maximum absolute residual change per time step.
// ============================================================================
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

// ============================================================================
//      This method returns the relative L2 norm of the residual change per time step.
// ============================================================================
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

// ============================================================================
//      This method returns the maximal velocity in the domain.
// ============================================================================
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

// ============================================================================
//      This method returns the maximal Reynolds number in the domain.
// ============================================================================
double GKSMesh::getMaxRe()
{
    return  this->getMaxVelocity()*this->param.L / this->fluidParam.nu;
}


// ============================================================================
//      This method returns the maximal Ma number in the domain.
// ============================================================================
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


// ============================================================================
//      This method returns wether the residual changes are below the defined
//      thresholds.
// ============================================================================
bool GKSMesh::isConverged(ConservedVariable residual)
{
    bool flag = true;

    flag = flag && ( residual.rho  < this->param.convergenceCriterium[0] );
    flag = flag && ( residual.rhoU < this->param.convergenceCriterium[1] );
    flag = flag && ( residual.rhoV < this->param.convergenceCriterium[2] );
    flag = flag && ( residual.rhoE < this->param.convergenceCriterium[3] );

    return flag;
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          File Output
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================


// ============================================================================
//      This method writes VTK and .dat files of cells and interfaces
// ============================================================================
void GKSMesh::writeOutputFiles()
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

    if ( param.csvOutput == true )
    {
        ostringstream filename;
        filename << "out/result_" << this->iter << ".dat";
        writeResultFields(filename.str());
    }

    if ( param.csvOutput == true )
    {
        ostringstream filename;
        filename << "out/resultBoundaryFlux_" << this->iter << ".dat";
        writeResultBoundaryFluxes(filename.str());
    }
}

// ============================================================================
//      This method writes a file with an overview of the simulation
//
//      Parameters:
//          filename:   filename
// ============================================================================
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

// ============================================================================
//      This method writes the time step size at the output intervals to a file.
//
//      Parameters:
//          filename:   filename
// ============================================================================
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

// ============================================================================
//      This method writes the times of the output intervals to a file.
//
//      Parameters:
//          filename:   filename
// ============================================================================
void GKSMesh::writeTime(string filename)
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

    for (vector<double>::iterator i = this->timeList.begin(); i != this->timeList.end(); ++i)
        file << (*i) << "\n";

    file.close();

    cout << "done!" << endl;
}

// ============================================================================
//      This method writes the time evolution of the residual changes to a file.
//
//      Parameters:
//          filename:   filename
// ============================================================================
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

// ============================================================================
//      This method writes the current field data per cell to a VTK file.
//
//      Parameters:
//          filename:   filename
// ============================================================================
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

// ============================================================================
//      This method writes the current field data per interface to a VTK file.
//
//      Parameters:
//          filename:   filename
// ============================================================================
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

// ============================================================================
//      This method writes the cell data of the current time step into a text
//      file for postprocessing for example in Matlab or Gnuplot. 
//
//      The order of the entries is:
//          1: CellCenter X
//          2: CellCenter Y
//          3: rho
//          4: U
//          5: V
//          6: T
//          7: Ghost Cell
//          8: rhoU
//          9: rhoV
//         10: rhoE
//         11: Cell Volume
//
//      Parameters:
//          filename:   filename
// ============================================================================
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
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getCenter().x << " ";
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getCenter().y << " ";
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getPrim().rho << " ";
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getPrim().U << " ";
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getPrim().V << " ";
        file << right << scientific << setw(25) << setfill(' ') << 1.0 / ( 2.0 * ( *i )->getPrim().L * this->fluidParam.R ) << " ";
        file << right <<               setw(3)  << setfill(' ') << ( *i )->isGhostCell() << " ";
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getCons().rhoU << " ";
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getCons().rhoV << " ";
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getCons().rhoE << " ";
        file << right << scientific << setw(25) << setfill(' ') << ( *i )->getVolume() << " ";
        file << endl;
    }

    file.close();

    cout << "done!" << endl;
}

// ============================================================================
//      This method writes the FluxDensities of the boundary interfaces into a
//      test file for postprocessing for example in Matlab or Gnuplot.
//
//      The order of the entries is:
//          1: InterfaceCenter X
//          2: InterfaceCenter Y
//          3: Normal X
//          4: Normal Y
//          5: F_rho
//          6: F_rhoU
//          7: F_rhoV
//          8: F_rhoE
//          9: F_sign
//         10: Interface Area
//
//      Parameters:
//          filename:   filename
// ============================================================================
void GKSMesh::writeResultBoundaryFluxes(string filename)
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

    for (vector<Interface*>::iterator i = this->InterfaceList.begin(); i != this->InterfaceList.end(); ++i)
    {
        if( (*i)-> isBoundaryInterface() )
        {
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getCenter().x << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getCenter().y << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getNormal().x << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getNormal().y << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getFluxDensity().rho  << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getFluxDensity().rhoU << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getFluxDensity().rhoV << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getFluxDensity().rhoE << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getFluxSign( (*i)->getCellInDomain() ) << " ";
            file << right << scientific << setw(25) << setfill(' ') << ( *i )->getArea() << " ";
            file << endl;
        }
    }

    file.close();

    cout << "done!" << endl;
}

// ============================================================================
//      This method export the mesh geometry into the Gambit neutral mesh
//      format for usage in other solvers.
//
//      Parameters:
//          filename:   filename
// ============================================================================
void GKSMesh::writeGambitNeutralFile(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
    // open file stream
    ofstream file;
    file.precision(11);
    file.open(filename.c_str());

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return;
    }
    
    // ========================================================================
    // count number of Nodes and Cells (without Ghost Cells)
    // ========================================================================
    unsigned long int numberOfNodes = 0;
    for(int i = 0; i < this->NodeList.size(); i++)  if(  this->NodeList[i]->ID != 9999999 ) numberOfNodes++;

    unsigned long int numberOfCells = 0;
    for(int i = 0; i < this->CellList.size(); i++)  if( !this->CellList[i]->isGhostCell() ) numberOfCells++;
    // ========================================================================
    
    // ========================================================================
    file << "        CONTROL INFO 2.4.6"                                    << endl;
    file << "** GAMBIT NEUTRAL FILE"                                        << endl;
    file << "SineDistortedMesh"                                             << endl;
    file << "PROGRAM:               GKSLenz     VERSION:  2.4.6"            << endl;
    file << " "                                                             << endl;
    file << "     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL"  << endl;
    file << right << setw(10) << setfill(' ') << numberOfNodes;
    file << right << setw(10) << setfill(' ') << numberOfCells;
    file << right << setw(10) << setfill(' ') << "1";
    file << right << setw(10) << setfill(' ') << "4";
    file << right << setw(10) << setfill(' ') << "2";
    file << right << setw(10) << setfill(' ') << "2";
    file << endl;
    file << "ENDOFSECTION"                                                  << endl;
    // ========================================================================

    // ========================================================================
    file << "   NODAL COORDINATES 2.4.6"                                    << endl;
    for(int i = 0; i < this->NodeList.size(); i++)
    {
        if(this->NodeList[i]->ID != 9999999)
        {
            file        << right <<               setw(10) << setfill(' ') << this->NodeList[i]->ID;
            file << " " << right << scientific << setw(19) << setfill(' ') << this->NodeList[i]->x;
            file << " " << right << scientific << setw(19) << setfill(' ') << this->NodeList[i]->y;
            file << endl;
        }
    }
    file << "ENDOFSECTION"                                                  << endl;
    // ========================================================================
    
    // ========================================================================
    file << "      ELEMENTS/CELLS 2.4.6"                                    << endl;
    for(int i = 0; i < this->CellList.size(); i++)
    {
        if( !this->CellList[i]->isGhostCell() )
        {
            file << right << setw(8) << setfill(' ') << this->CellList[i]->getID();
            file << right << setw(3) << setfill(' ') << "2";
            file << right << setw(3) << setfill(' ') << "4";
            file << right << setw(1) << setfill(' ') << " ";
            file << right << setw(8) << setfill(' ') << this->CellList[i]->getNode(0).ID;
            file << right << setw(8) << setfill(' ') << this->CellList[i]->getNode(1).ID;
            file << right << setw(8) << setfill(' ') << this->CellList[i]->getNode(2).ID;
            file << right << setw(8) << setfill(' ') << this->CellList[i]->getNode(3).ID;
            file << endl;
        }
    }
    file << "ENDOFSECTION"                                                  << endl;
    // ========================================================================
    
    // ========================================================================
    file << "       ELEMENT GROUP 2.4.6"                                    << endl;
    file << "GROUP:          1 ELEMENTS: ";
    file << right << setw(10) << setfill(' ') << numberOfCells;
    file << " MATERIAL:          2 NFLAGS:          1"                      << endl;
    file << "                           fluid"                              << endl;
    file << "       0"                                                      << endl;
    int j = 0;
    for(int i = 0; i < this->CellList.size(); i++)
    {
        if( !this->CellList[i]->isGhostCell() )
        {
            file << right << setw(8) << setfill(' ') << this->CellList[i]->getID();
            j++;
            if(j%10 == 0) file << endl;
        }
    }
    file << endl;
    file << "ENDOFSECTION"                                                  << endl;
    // ========================================================================
    
    // ================================================================================================================
    // ================================================================================================================
    // ================================================================================================================
    // ================================================================================================================
    
    // ========================================================================
    int numberBC = 0;
    for(int i = 0; i < this->InterfaceList.size(); i++)
        if( this->InterfaceList[i]->getBoundaryCondition() == this->BoundaryConditionList[0] )  numberBC++;

    file << " BOUNDARY CONDITIONS 2.4.6"                                    << endl;
    file << "                           left       1";
    file << right << setw(8) << setfill(' ') << numberBC;
    file << "       0       6" << endl;
    for(int i = 0; i < this->InterfaceList.size(); i++)
    {
        if( this->InterfaceList[i]->getBoundaryCondition() == this->BoundaryConditionList[0] )
        {
            file << right << setw(10) << setfill(' ') << this->InterfaceList[i]->getCellInDomain()->getID();
            file << right << setw(5)  << setfill(' ') << 2;
            file << right << setw(5)  << setfill(' ') << 4;
            file << endl;
        }
    }
    file << "ENDOFSECTION"                                                  << endl;
    // ========================================================================
    
    // ========================================================================
    numberBC = 0;
    for(int i = 0; i < this->InterfaceList.size(); i++)
        if( this->InterfaceList[i]->getBoundaryCondition() == this->BoundaryConditionList[1] )  numberBC++;

    file << " BOUNDARY CONDITIONS 2.4.6"                                    << endl;
    file << "                         bottom       1";
    file << right << setw(8) << setfill(' ') << numberBC;
    file << "       0       6" << endl;
    for(int i = 0; i < this->InterfaceList.size(); i++)
    {
        if( this->InterfaceList[i]->getBoundaryCondition() == this->BoundaryConditionList[1] )
        {
            file << right << setw(10) << setfill(' ') << this->InterfaceList[i]->getCellInDomain()->getID();
            file << right << setw(5)  << setfill(' ') << 2;
            file << right << setw(5)  << setfill(' ') << 1;
            file << endl;
        }
    }
    file << "ENDOFSECTION"                                                  << endl;
    // ========================================================================
    
    // ========================================================================
    numberBC = 0;
    for(int i = 0; i < this->InterfaceList.size(); i++)
        if( this->InterfaceList[i]->getBoundaryCondition() == this->BoundaryConditionList[2] )  numberBC++;

    file << " BOUNDARY CONDITIONS 2.4.6"                                    << endl;
    file << "                          right       1";
    file << right << setw(8) << setfill(' ') << numberBC;
    file << "       0       6" << endl;
    for(int i = 0; i < this->InterfaceList.size(); i++)
    {
        if( this->InterfaceList[i]->getBoundaryCondition() == this->BoundaryConditionList[2] )
        {
            file << right << setw(10) << setfill(' ') << this->InterfaceList[i]->getCellInDomain()->getID();
            file << right << setw(5)  << setfill(' ') << 2;
            file << right << setw(5)  << setfill(' ') << 2;
            file << endl;
        }
    }
    file << "ENDOFSECTION"                                                  << endl;
    // ========================================================================
    
    // ========================================================================
    numberBC = 0;
    for(int i = 0; i < this->InterfaceList.size(); i++)
        if( this->InterfaceList[i]->getBoundaryCondition() == this->BoundaryConditionList[3] )  numberBC++;

    file << " BOUNDARY CONDITIONS 2.4.6"                                    << endl;
    file << "                            top       1";
    file << right << setw(8) << setfill(' ') << numberBC;
    file << "       0       6" << endl;
    for(int i = 0; i < this->InterfaceList.size(); i++)
    {
        if( this->InterfaceList[i]->getBoundaryCondition() == this->BoundaryConditionList[3] )
        {
            file << right << setw(10) << setfill(' ') << this->InterfaceList[i]->getCellInDomain()->getID();
            file << right << setw(5)  << setfill(' ') << 2;
            file << right << setw(5)  << setfill(' ') << 3;
            file << endl;
        }
    }
    file << "ENDOFSECTION"                                                  << endl;
    // ========================================================================

}

// ============================================================================
//      This private method writes the Cell geometry of the mesh into the 
//      specified File in the VTK format including the file header.
//
//      Parameters:
//          file:   destination file
// ============================================================================
void GKSMesh::writeCellGeometry(ofstream& file)
{
    
    int numberOfCells = 0;
    for(int i = 0; i < this->CellList.size(); i++)  
        if( !this->CellList[i]->isGhostCell() || this->param.ghostOutput)  numberOfCells++;

    // write VTK Header
    file << "# vtk DataFile Version 1.0\n";
    file << "by Stephan Lenz\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // write nodes
    //( one dummy node with the ID 0 must be written )
    file << "POINTS " << 4 * numberOfCells + 1 << " double\n";
    file << "0.0 0.0 0.0 \n";

    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << (*i)->writeNodes();
    }

    // write elements
    file << "CELLS " << numberOfCells << " " << 5 * numberOfCells << endl;
    for (int i = 0; i < this->CellList.size(); ++i)
    {
        if( !CellList[i]->isGhostCell() || this->param.ghostOutput )
            file << 4 << " " << i * 4 + 1
                << " " << i * 4 + 2
                << " " << i * 4 + 3
                << " " << i * 4 + 4 << endl;
    }

    // write element tyes( 9 = quad element )
    file << "CELL_TYPES " << numberOfCells << endl;
    for (int i = 0; i < this->CellList.size(); ++i)
    {
        if( !CellList[i]->isGhostCell() || this->param.ghostOutput )
            file << "9" << endl;
    }
}

// ============================================================================
//      This private method writes the Interface geometry of the mesh into the 
//      specified File in the VTK format including the file header.
//
//      Parameters:
//          file:   destination file
// ============================================================================
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

// ============================================================================
//      This private method writes the Cell data into the 
//      specified File in the VTK format.
//
//      Parameters:
//          file:   destination file
// ============================================================================
void GKSMesh::writeCellData(ofstream& file)
{
    int numberOfFields = 14;
    if ( this->param.resOutput )
        numberOfFields += 4;
    
    int numberOfCells = 0;
    for(int i = 0; i < this->CellList.size(); i++)  
        if( !this->CellList[i]->isGhostCell() || this->param.ghostOutput)  numberOfCells++;

    // write cell data ( ID and stress )
    file << "CELL_DATA " << numberOfCells << endl;
    file << "FIELD Lable " << numberOfFields << "\n";

    // ================================================================================================================

    file << "rho 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getPrim().rho << endl;
    }

    file << "U 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getPrim().U << endl;
    }

    file << "V 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getPrim().V << endl;
    }

    file << "Lambda 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getPrim().L << endl;
    }

    // ================================================================================================================

    file << "GhostCell 1 " << numberOfCells << " int\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
        {
            if ( ( *i )->isGhostCell() )
                file << 1 << endl;
            else
                file << 0 << endl;
        }
    }

    // ================================================================================================================

    file << "rhoU 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getCons().rhoU << endl;
    }

    file << "rhoV 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getCons().rhoV << endl;
    }

    file << "rhoE 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getCons().rhoE << endl;
    }

    // ================================================================================================================

    file << "p 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getPrim().rho / ( 2.0 * ( *i )->getPrim().L ) << endl;
    }

    file << "T 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << 1.0 / ( 2.0 * this->fluidParam.R * ( *i )->getPrim().L ) << endl;
    }

    // ================================================================================================================

    if ( this->param.resOutput )
    {
        file << "res_rho 1 " << numberOfCells << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            if( !(*i)->isGhostCell() || this->param.ghostOutput )
                file << ( *i )->getLocalResidual().rho << endl;
        }

        file << "res_rhoU 1 " << numberOfCells << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            if( !(*i)->isGhostCell() || this->param.ghostOutput )
                file << ( *i )->getLocalResidual().rhoU << endl;
        }

        file << "res_rhoV 1 " << numberOfCells << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            if( !(*i)->isGhostCell() || this->param.ghostOutput )
                file << ( *i )->getLocalResidual().rhoV << endl;
        }

        file << "res_rhoE 1 " << numberOfCells << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            if( !(*i)->isGhostCell() || this->param.ghostOutput )
                file << ( *i )->getLocalResidual().rhoE << endl;
        }
    }
    // ================================================================================================================
    file << "update_rho 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getUpdate().rho << endl;
    }

    file << "update_rhoU 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getUpdate().rhoU << endl;
    }

    file << "update_rhoV 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getUpdate().rhoV << endl;
    }

    file << "update_rhoE 1 " << numberOfCells << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << ( *i )->getUpdate().rhoE << endl;
    }
    // ================================================================================================================
    for(int j = 0; j < 4; j++)
    {
        file << "VECTORS Connectivity" << j << " double\n";
        for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
        {
            if( !(*i)->isGhostCell() || this->param.ghostOutput )
                file << (*i)->getConnectivity(j).x << " " << (*i)->getConnectivity(j).y << " 0.0" << endl;
        }
    }
    // ================================================================================================================
    file << "VECTORS dRho" << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << (*i)->getGradientX().rho << " " << (*i)->getGradientY().rho << " 0.0" << endl;
    }
    file << "VECTORS dRhoU" << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << (*i)->getGradientX().rhoU << " " << (*i)->getGradientY().rhoU << " 0.0" << endl;
    }
    file << "VECTORS dRhoV" << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << (*i)->getGradientX().rhoV << " " << (*i)->getGradientY().rhoV << " 0.0" << endl;
    }
    file << "VECTORS dRhoE" << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        if( !(*i)->isGhostCell() || this->param.ghostOutput )
            file << (*i)->getGradientX().rhoE << " " << (*i)->getGradientY().rhoE << " 0.0" << endl;
    }
    // ================================================================================================================
}

// ============================================================================
//      This private method writes the Interface data into the 
//      specified File in the VTK format.
//
//      Parameters:
//          file:   destination file
// ============================================================================
void GKSMesh::writeInterfaceData(ofstream & file)
{
    // write cell data ( ID and stress )
    file << "POINT_DATA " << this->InterfaceList.size() << endl;
    file << "FIELD Lable 8\n";

    file << "F_rho 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << ( *i )->getFluxDensity().rho << endl;
    }

    file << "F_rhoU 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << ( *i )->getFluxDensity().rhoU << endl;
    }

    file << "F_rhoV 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        file << ( *i )->getFluxDensity().rhoV << endl;
    }

    file << "F_rhoE 1 " << this->InterfaceList.size() << " double\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
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

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================