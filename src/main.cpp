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

#include "GKSSolverPush.h"
#include "BoundaryCondition.h"
#include "mshReader.h"
#include "outputWriter.h"
#include <iostream>
#include <sstream>
#include <chrono>

# define M_PI           3.14159265358979323846  /* pi */

using namespace std;

int main(int argc, char* argv[])
{
    // ================================================================================================================
    // ================================================================================================================
    // ================================================================================================================
    // ================================================================================================================
    // ================================================================================================================

    // These loops can be used for Convergence studies
    double ReList[] = {40.0, 100.0, 400.0, 1000.0};
    int    nyList[] = {8, 16, 32, 64, 128};
    double RaList[] = {1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8};
    //for (int i = 0; i < 4;i++)      //ReList
    //for (int j = 0; j < 5;j++)      // nyList
    //for(int i = 0; i < 6; i++)      // RaList
    if(true){
        ///*
        
        // ============================================================================================================
        // ============================================================================================================
        // ============================================================================================================
        //
        //              Problem-Definition:
        //                  Poiseuille-Flow (Force driven)
        //
        // ============================================================================================================
        // ============================================================================================================
        // ============================================================================================================

        // ========================================================================
        //                  Simulation parameters
        // ========================================================================

        Parameters param;

        param.verbose     = false;                      // detailed screen output
        param.fluxOutput  = false;                      // VTK files for interfaces
        param.resOutput   = false;                      // include residuals in VTK files
        param.ghostOutput = true;                      // include ghost cells in VTK files
        param.csvOutput   = false;                       // output csv files for postprocessing

        param.numberOfIterations = 200000000;            // maximal number of Iterations
        param.outputIntervalVTK  = 10000;              // Output interval for VTK Files (and .dat files)
        param.outputInterval     = 1000;              // Output interval for Output on the screen

        // Abortion criteria for the different conserved variables
        // These are thresholds for relative residual changes
        // ||W^n+1 - W^n||_L2 / ||W^n||_L2
        param.convergenceCriterium[0] = 1.0;
        param.convergenceCriterium[1] = 1.0e-10;
        param.convergenceCriterium[2] = 1.0;
        param.convergenceCriterium[3] = 1.0;

        param.CFL = 0.01;                                // CFL number for time step computation
        
        // ========================================================================
        //                  Fluid and domain parameters
        // ========================================================================

        FluidParameter fluidParam;
        param.L = 1;          // reference length for Re number

        fluidParam.K = 1.0;                                                   
        fluidParam.nu = 0.001;//(u0*param.L)/Re;                              
        fluidParam.R = 200.0;
        fluidParam.Pr = 1.0;

        fluidParam.rhoReference = 1.0;
        fluidParam.uReference   = 1.5;
        fluidParam.vReference   = 0.0;
        fluidParam.lambdaReference = 1.0 / (2.0 * fluidParam.R * 1.0);

        fluidParam.Force.x = 0.1;//(u0*8.0*fluidParam.nu) / (param.L*param.L);
        fluidParam.Force.y = 0.0;
        fluidParam.BoussinesqForce.x = 0.0;
        fluidParam.BoussinesqForce.y = 0.0;
        
        // ============================================================================================================
        // ============================================================================================================
        // ============================================================================================================
        //
        //              Start the simulation and output several files
        //
        // ============================================================================================================
        // ============================================================================================================
        // ============================================================================================================


        GKSSolver* solverPush = new GKSSolverPush(param, fluidParam);

        if( ! solverPush->readProblem("msh/TurekBenchmark.mesh.fine") )
        //if( ! solverPush->readProblem("msh/TurekPeriodic/TurekPeriodic.mesh") )
        //if( ! solverPush->readProblem("msh/Channel/Channel_64x32") )
        {
            system("pause");
            return false;
        }

        outputWriter::writeCellVTK("out/Mesh.vtk", *solverPush);
        outputWriter::writeFaceVTK("out/Connectivity.vtk", *solverPush);

        solverPush->iterate();

        system("pause");
        delete solverPush;
    }
}