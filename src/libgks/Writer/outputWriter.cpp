#include "libgks/Writer/outputWriter.h"
#include <iostream>

void outputWriter::writeCellVTK(string filename, GKSSolver& solver)
{
	ofstream file;
	open(file, filename + string(".vtk"));
    
    // write VTK Header
    file << "# vtk DataFile Version 1.0\n";
    file << "by Stephan Lenz\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    
    file << "POINTS " << solver.getNumberOfNodes() << " double\n";
    for(int node = 0; node < solver.getNumberOfNodes(); ++node)
    {
        file << solver.getNode(node).x << " " << solver.getNode(node).y << " 0.0" << endl;
    }

    file << "CELLS " << solver.getNumberOfCells() << " " << 5 * solver.getNumberOfCells() << endl;
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << "4 " << solver.getCell2Node(cell)[0] << " "
                     << solver.getCell2Node(cell)[1] << " "
                     << solver.getCell2Node(cell)[2] << " ";

        if( -1 == solver.getCell2Node(cell)[3] )
            file << solver.getCell2Node(cell)[2] << endl;
        else
            file << solver.getCell2Node(cell)[3] << endl;
    }

    file << "CELL_TYPES " << solver.getNumberOfCells() << endl;
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << "9" << endl;
    }
    
    file << "CELL_DATA " << solver.getNumberOfCells() << endl;
    file << "FIELD Lable " << 8 << "\n";
    file << "CellID 1 " << solver.getNumberOfCells() << " int\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << cell << endl;
    }
    file << "rho 1 " << solver.getNumberOfCells() << " double\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << solver.getPrim(cell).rho << endl;
    }
    file << "U 1 " << solver.getNumberOfCells() << " double\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << solver.getPrim(cell).U << endl;
    }
    file << "V 1 " << solver.getNumberOfCells() << " double\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << solver.getPrim(cell).V << endl;
    }
    file << "T 1 " << solver.getNumberOfCells() << " double\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << 1.0 / ( 2.0 * solver.getFluidParam().R * solver.getPrim(cell).L ) << endl;
    }
    file << "p 1 " << solver.getNumberOfCells() << " double\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << solver.getPrim(cell).rho / ( 2.0 * solver.getPrim(cell).L ) << endl;
    }
    file << "Lambda 1 " << solver.getNumberOfCells() << " double\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << solver.getPrim(cell).L << endl;
    }
    file << "BC 1 " << solver.getNumberOfCells() << " int\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << solver.isGhostCell(cell) << endl;
    }

    file << "VECTORS Velocity double\n";
    for (int cell = 0; cell < solver.getNumberOfCells(); ++cell)
    {
        file << solver.getPrim( cell ).U << " " << solver.getPrim( cell ).V << " 0.0" << endl;
    }

    file.close();

    cout << "done!" << endl;
}

void outputWriter::writeFaceVTK(string filename, GKSSolver& solver)
{
	ofstream file;
    open( file, filename );
    
    // write VTK Header
    file << "# vtk DataFile Version 1.0\n";
    file << "by Stephan Lenz\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    
    file << "POINTS " << solver.getNumberOfInterfaces() << " double\n";
    for(int face = 0; face < solver.getNumberOfInterfaces(); ++face)
    {
        file << solver.getInterfaceCenter(face).x << " " << solver.getInterfaceCenter(face).y << " 0.0" << endl;
    }

    file << "CELLS " << solver.getNumberOfInterfaces() << " " << 2 * solver.getNumberOfInterfaces() << endl;
    for (int cell = 0; cell < solver.getNumberOfInterfaces(); ++cell)
    {
        file << "1 " << cell << endl;
    }

    file << "CELL_TYPES " << solver.getNumberOfInterfaces() << endl;
    for (int cell = 0; cell < solver.getNumberOfInterfaces(); ++cell)
    {
        file << "1" << endl;
    }

    file << "POINT_DATA " << solver.getNumberOfInterfaces() << endl;
    file << "FIELD Lable 2\n";

    file << "ID 1 " << solver.getNumberOfInterfaces() << " double\n";
    for ( int face = 0; face < solver.getNumberOfInterfaces(); ++face )
    {
        file << face << endl;
    }

    file << "Area 1 " << solver.getNumberOfInterfaces() << " double\n";
    for ( int face = 0; face < solver.getNumberOfInterfaces(); ++face )
    {
        file << solver.getInterfaceArea(face) << endl;
    }

    file << "VECTORS normal double\n";
    for ( int face = 0; face < solver.getNumberOfInterfaces(); ++face )
    {
        file << solver.getInterfaceNormal(face).x << " " << solver.getInterfaceNormal(face).y << " 0.0" << endl;
    }

    file << "VECTORS posCell double\n";
    for ( int face = 0; face < solver.getNumberOfInterfaces(); ++face )
    {
        file << solver.getCellCenter( solver.getPosCell(face) ).x - solver.getInterfaceCenter(face).x << " " 
             << solver.getCellCenter( solver.getPosCell(face) ).y - solver.getInterfaceCenter(face).y << " 0.0" << endl;
    }

    file << "VECTORS negCell double\n";
    for ( int face = 0; face < solver.getNumberOfInterfaces(); ++face )
    {
        file << solver.getCellCenter( solver.getNegCell(face) ).x - solver.getInterfaceCenter(face).x << " " 
             << solver.getCellCenter( solver.getNegCell(face) ).y - solver.getInterfaceCenter(face).y << " 0.0" << endl;
    }

    file.close();

    cout << "done!" << endl;
}

void outputWriter::initFile(string filename)
{
    ofstream file;
    file.open( filename.c_str(), fstream::trunc ); 
    file.close();
}

bool outputWriter::open(ofstream& file, string filename)
{
    std::cout << "Wrinting file " << filename << " ... ";
	
	file.open(filename.c_str());
    file.precision(15);

	if (!file.is_open()) {
		std::cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return false;
	}

    return true;
}

void outputWriter::writeOverview(string filename, GKSSolver & solver)
{
    ofstream file;
    open( file, filename + string(".Overview.dat" ) );

    file << " ========== Fluid Parameters ==========";
    file << endl;
    file << "nu =\t " << solver.getFluidParam().nu << endl;
    file << "K  =\t " << solver.getFluidParam().K << endl;
    file << "R  =\t " << solver.getFluidParam().R << endl;
    file << "Fx =\t " << solver.getFluidParam().Force.x << endl;
    file << "Fy =\t " << solver.getFluidParam().Force.y << endl;
    file << "Pr =\t " << solver.getFluidParam().Pr << endl;
    file << endl;

    file << " ========== Simulation Parameters ==========";
    file << endl;
    file << "Max Number of Iteratios:               " << solver.getParameters().numberOfIterations << endl;
    file << "VTK-File Output Interval:              " << solver.getParameters().outputIntervalVTK << endl;
    file << "Convergence History Output Interval:   " << solver.getParameters().outputInterval << endl;
    file << "Maximal Simulation Time:               " << solver.getParameters().maxTime << endl;
    file << "Convergence Criterium:               ( " << solver.getParameters().convergenceCriterium.rho  << ", "
                                                      << solver.getParameters().convergenceCriterium.rhoU << ", "
                                                      << solver.getParameters().convergenceCriterium.rhoV << ", "
                                                      << solver.getParameters().convergenceCriterium.rhoE << " )" << endl;
    file << "CFL =\t" << solver.getParameters().CFL << endl;
    file << endl;

    file << " ========== Flow Characteristics ==========";
    file << endl;
    file << "Umax = " << solver.getMaxVelocity() << " m/s" << endl;
    file << "Ma   = " << solver.getMaxMa() << endl;
    file << endl;

    file << " ========== Simulation Results ==========";
    file << endl;
    if ( solver.isConverged( solver.getL2GlobalResidual() ) )
    {
        file << "Simulation Converged!" << endl;
        file << "Time steps: " << solver.getIter() << endl;
    }
    else
    {
        file << "Simulation did not converge" << endl;
    }
    file << "Final Residuals:" << endl;
    ConservedVariable residual = solver.getL2GlobalResidual();
    file << "r_rho  = " << residual.rho << "\t ";
    file << "r_rhoU = " << residual.rhoU << "\t ";
    file << "r_rhoV = " << residual.rhoV << "\t ";
    file << "r_rhoE = " << residual.rhoE << "\t ";
    file << endl;
    file << endl;
    file << "Real Time simulated : " << solver.getTime() << " s" << endl;
    file << "Time to solution:     " << solver.getComputationTime() << " s" << endl;
    file << "CellUpdatesPerSecond: " << solver.getCellUpdatesPerSecond() << endl;
    file.close();

cout << "done!" << endl;
}

void outputWriter::writeConvergenceHistory(string filename, GKSSolver& solver)
{
    ofstream file;
    file.open( ( filename + string(".ConvHist.dat" ) ).c_str(), std::fstream::app );

    ConservedVariable residual = solver.getL2GlobalResidual();
    file << residual.rho  << ' '
         << residual.rhoU << ' '
         << residual.rhoV << ' '
         << residual.rhoE << endl;

    file.close();
}
