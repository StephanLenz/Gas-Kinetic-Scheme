
#include "GKSMesh.h"
#include "Cell.h"
#include "Types.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>    //min()
#include <chrono>

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

void GKSMesh::generateRectMesh(double lengthX, double lengthY, int nx, int ny)
{
	double dx = lengthX / (double)nx;
	double dy = lengthY / (double)ny;

	this->lengthX = lengthX;
	this->lengthY = lengthY;

	Cell*		tmpCell;
	Interface*  tmpInterface;
    float2      normal;
    float2      center;

	//=========================================================================
	//=========================================================================
	//		Cell generation
	//			including ghost cells
	//=========================================================================
	//=========================================================================
    BoundaryCondition* currentBC = NULL;
	for (int i = -1; i < ny + 1; i++)       // Y-Direction
	{

		for (int j = -1; j < nx + 1; j++)   // X-Direction
		{
            if (i == -1)         currentBC = BoundaryConditionList[1];
            else if (i == ny)    currentBC = BoundaryConditionList[3];
            else if (j == -1)    currentBC = BoundaryConditionList[0];
            else if (j == nx)    currentBC = BoundaryConditionList[2];
            else                 currentBC = NULL;

			//                      cell centerX         cell centerY
			tmpCell = new Cell(((double)j + 0.5)*dx, ((double)i + 0.5)*dy, dx, dy, currentBC, this->fluidParam);
			// add interface to list
			this->CellList.push_back(tmpCell);
		}
	}

	//=========================================================================
	//=========================================================================
	//						F interface generation
	//=========================================================================
	//=========================================================================
    normal.x = 1;
    normal.y = 0;
	for (int i = 0; i <= ny + 1; i++)       // Y-Direction
	{
		for (int j = 0; j < nx + 1; j++)    // X-Direction
		{
            center.x = (double)j*dx;
            center.y = ( (double)i - 0.5 )*dy;

            Cell* posCell = this->CellList[i*(nx + 2) + (j + 1)];
            Cell* negCell = this->CellList[i*(nx + 2) + j];

			// create a new interface with the adjacent cells
			tmpInterface = new Interface(negCell, posCell, center, normal, this->fluidParam);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}

	//=========================================================================
	//=========================================================================
	//						G interface generation
	//=========================================================================
	//=========================================================================
    normal.x = 0;
    normal.y = 1;
	for (int i = 0; i < ny + 1; i++)        // Y-Direction
	{
		for (int j = 0; j <= nx + 1; j++)   // X-Direction
		{
            center.x = ( (double)j - 0.5 )*dx;
            center.y = (double)i*dy;

            Cell* posCell = this->CellList[(i + 1)*(nx + 2) + j];
            Cell* negCell = this->CellList[i*(nx + 2) + j];

			// create a new interface with the adjacent cells
			tmpInterface = new Interface(negCell, posCell, center, normal, this->fluidParam);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}


	return;
}

void GKSMesh::generateRectMeshPeriodic(double lengthX, double lengthY, int nx, int ny)
{
    double dx = lengthX / (double)nx;
    double dy = lengthY / (double)ny;

    this->lengthX = lengthX;
    this->lengthY = lengthY;

    Cell*		tmpCell;
    Interface*  tmpInterface;
    float2      normal;
    float2      center;

    //=========================================================================
    //=========================================================================
    //		Cell generation
    //			including ghost cells in y-direction
    //=========================================================================
    //=========================================================================
    BoundaryCondition* currentBC = NULL;
    for (int i = -1; i < ny + 1; i++)       // Y-Direction
    {

        for (int j = 0; j < nx; j++)   // X-Direction
        {
            if (i == -1)         currentBC = BoundaryConditionList[0];
            else if (i == ny)    currentBC = BoundaryConditionList[1];
            else                 currentBC = NULL;

            //                      cell centerX         cell centerY
            tmpCell = new Cell(((double)j + 0.5)*dx, ((double)i + 0.5)*dy, dx, dy, currentBC, this->fluidParam);
            // add interface to list
            this->CellList.push_back(tmpCell);
        }
    }

    //=========================================================================
    //=========================================================================
    //						F interface generation
    //=========================================================================
    //=========================================================================
    normal.x = 1;
    normal.y = 0;
    for (int i = 0; i <= ny + 1; i++)       // Y-Direction
    {
        for (int j = 0; j < nx; j++)    // X-Direction
        {
            center.x = (double)j * dx;
            center.y = ( (double)i - 0.5 )*dy;

            Cell* negCell;
            Cell* posCell;

            if (j == 0)
                negCell = CellList[i*(nx) + (nx-1)];
            else
                negCell = CellList[i*(nx) + (j-1)];

            if (j == nx)
                posCell = CellList[i*(nx)];
            else
                posCell = CellList[i*(nx) + j];

            // create a new interface with the adjacent cells
            tmpInterface = new Interface( negCell, posCell, center, normal, this->fluidParam);
            // add itnerface to list
            this->InterfaceList.push_back(tmpInterface);
        }
    }

    //=========================================================================
    //=========================================================================
    //						G interface generation
    //=========================================================================
    //=========================================================================
    normal.x = 0;
    normal.y = 1;
    for (int i = 0; i < ny + 1; i++)        // Y-Direction
    {
        for (int j = 0; j < nx; j++)        // X-Direction
        {
            center.x = ( (double)j + 0.5 ) * dx;
            center.y = (double)i * dy;

            Cell* posCell = CellList[( i + 1 )*(nx)+j];
            Cell* negCell = CellList[i*(nx)+j];

            // create a new interface with the adjacent cells
            tmpInterface = new Interface(negCell, posCell, center, normal, this->fluidParam);
            // add itnerface to list
            this->InterfaceList.push_back(tmpInterface);
        }
    }

    return;
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

void GKSMesh::addBoundaryCondition( int rhoType, int UType, int VType, int TType, 
                                    double rho, double U, double V, double T)
{
    BoundaryCondition* tmp = new BoundaryCondition( rhoType, UType, VType, TType,
                                                    rho, U, V, T);
    BoundaryConditionList.push_back(tmp);
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

ConservedVariable GKSMesh::getGlobalResidual()
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

void GKSMesh::timeStep()
{
    this->iter++;

    if (this->param.verbose) cout << "Iterration: " << this->iter << endl;

    if(this->param.verbose) cout << "  Compute Timestep ..." << endl;
    this->computeGlobalTimestep();
    if (this->param.verbose) cout << "    dt = " << this->dt << endl;

    this->dtList.push_back(this->dt);

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
            ConservedVariable residual = this->getGlobalResidual();

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
        //cout << this->toString();
        //cout << this->cellValuesToString();
    }
    // ========================================================================
    // ========================================================================
    // ========================================================================
    
    chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();
    this->computationTime = chrono::duration_cast<chrono::seconds>( endTime - startTime ).count();

    cout << "Time to Solution: " << this->computationTime << " s" << endl;
}

bool GKSMesh::isConverged(ConservedVariable residual)
{
    bool flag = true;

    flag = flag && ( residual.rho  < this->param.convergenceCriterium );
    flag = flag && ( residual.rhoU < this->param.convergenceCriterium );
    flag = flag && ( residual.rhoV < this->param.convergenceCriterium );
    flag = flag && ( residual.rhoE < this->param.convergenceCriterium );

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
    file << "Convergence Criterium:               " << this->param.convergenceCriterium << endl;
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
            // check wether the profile location x is located in this cell
            if ( fabs( ( *i )->getCenter().x - x ) <= 0.5 * ( *i )->getDx().x )
            {
                file << ( *i )->getCenter().y << " " << ( *i )->getPrim().U << "\n";
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
    file << "POINTS " << 4 * this->CellList.size() + 1 << " float\n";
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
    file << "POINTS " << this->InterfaceList.size()<< " float\n";

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
    // write cell data ( ID and stress )
    file << "CELL_DATA " << this->CellList.size() << endl;
    file << "FIELD Lable 13\n";

    // ================================================================================================================

    file << "rho 1 " << this->CellList.size() << " double\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getPrim().rho << endl;
    }

    file << "U 1 " << this->CellList.size() << " double\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getPrim().U << endl;
    }

    file << "V 1 " << this->CellList.size() << " double\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getPrim().V << endl;
    }

    file << "Lambda 1 " << this->CellList.size() << " double\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getPrim().L << endl;
    }

    // ================================================================================================================

    file << "GhostCell 1 " << this->CellList.size() << " int\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        if ((*i)->isGhostCell())
            file << 1 << endl;
        else
            file << 0 << endl;
    }

    // ================================================================================================================

    file << "rhoU 1 " << this->CellList.size() << " double\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getCons().rhoU << endl;
    }

    file << "rhoV 1 " << this->CellList.size() << " double\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getCons().rhoV << endl;
    }

    file << "rhoE 1 " << this->CellList.size() << " double\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getCons().rhoE << endl;
    }

    // ================================================================================================================

    file << "p 1 " << this->CellList.size() << " double\n";
    for ( vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i )
    {
        file << ( *i )->getPrim().rho / ( 2.0 * ( *i )->getPrim().L ) << endl;
    }

    // ================================================================================================================

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

    // ================================================================================================================
}

void GKSMesh::writeInterfaceData(ofstream & file)
{
    // write cell data ( ID and stress )
    file << "POINT_DATA " << this->InterfaceList.size() << endl;
    file << "FIELD Lable 5\n";

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

    file << "GhostInterface 1 " << this->InterfaceList.size() << " int\n";
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        if ((*i)->isGhostInterface())
            file << 1 << endl;
        else
            file << 0 << endl;
    }
}

