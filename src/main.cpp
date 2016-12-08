// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      main.cpp
//
//      This file defines the Problem setup and starts the simulation.
//
// ============================================================================

#include "outputWriter.h"
#include "GKSSolver.h"
#include "GKSSolverPush.h"
#include "paramReader.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    Parameters     param;
    FluidParameter fluidParam;

    if( argc == 1 )
    {
        cout << "Error: No Simulation Name given." << endl;
        system("pause");
        return 1;
    }

    param.simulationName = argv[1];

    if( ! paramReader::read( string("msh/") + param.simulationName + string(".gksparam"), param, fluidParam ) )
    {
        system("pause");
        return 1;
    }

    GKSSolver* solverPush = new GKSSolverPush(param, fluidParam);

    if( ! solverPush->readProblem( string("msh/") + param.simulationName ) )
    {
        system("pause");
        return 1;
    }

    outputWriter::writeCellVTK(string("out/") + param.simulationName + string(".mesh.vtk" ), *solverPush);
    outputWriter::writeFaceVTK(string("out/") + param.simulationName + string(".connectivity.vtk" ), *solverPush);

    solverPush->iterate();

    system("pause");
    delete solverPush;
    return 0;
}
