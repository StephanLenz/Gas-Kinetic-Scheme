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
        param.outputInterval     = 10000;              // Output interval for Output on the screen

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
        double Re = 10.0;        // Reynolds number
        double u0 = 0.1;        // Velocits in the mid of the channel
        double T  = 1.0;
        param.L = 1;          // reference length for Re number

        fluidParam.K = 1.0;                                                   // internal degrees of freedom
        fluidParam.nu = 0.001;//(u0*param.L)/Re;                                    // kinematic viskosity
        fluidParam.R = 200.0;                                               // specific gas constant
        fluidParam.Force.x = 0.0;//(u0*8.0*fluidParam.nu) / (param.L*param.L);    // acceleration in x direction [m/s^2]
        fluidParam.Force.y = 0.0;                                           // acceleration in y direction [m/s^2]
        fluidParam.BoussinesqForce.x = 0.0;                                 // acceleration only allpied to density variations [m/s^2]
        fluidParam.BoussinesqForce.y = 0.0;                                 // acceleration only allpied to density variations [m/s^2]
        fluidParam.rhoReference = 1.0;                                      // reference density
        fluidParam.uReference   = 1.5;
        fluidParam.vReference   = 0.0;
        fluidParam.lambdaReference = 1.0 / (2.0 * fluidParam.R * T);        // reference temperature
        fluidParam.Pr = 1.0;                                                // Prandl number 
        
        // ============================================================================================================
        // ============================================================================================================
        // ============================================================================================================
        //
        //              Start the simulation and output several files
        //
        // ============================================================================================================
        // ============================================================================================================
        // ============================================================================================================

        //cout << mesh->toString();

        //mesh->writeMeshAsText("out/Mesh.txt");

        //mesh->writeVTKFile("out/InitialState.vtk");
        //mesh->writeVTKFileFlux("out/InitialStateFlux.vtk");
        //mesh->writeGambitNeutralFile("out/SineDistortedMesh.neu");

        // ====================================================================
        // ====================================================================
        //              Run the actual simulation
        // ====================================================================
        // ====================================================================

        GKSSolver* solverPush = new GKSSolverPush(param, fluidParam);

        //if( ! solverPush->readMeshFromMshFile("msh/SquareQuadGradedPeriodicGhostIsothermalWall32.msh") )
        if( ! solverPush->readMeshFromMshFile("msh/TurekBenchmark/TurekBenchmark2G.mesh.3.msh") )
        {
            system("pause");
            return false;
        }

        solverPush->writeVTK("out/Mesh.vtk");
        solverPush->writeInterfaceVTK("out/Connectivity.vtk");

        //mesh->iterate();

        //solverPull->iterate();
        solverPush->iterate();
        //solverSOA->iterate();
        //solverAOS->iterate();

        //solverPull->writeDataToMeshObject(*mesh);
        //mesh->writeVTKFile( "out/solverPull.vtk" );

        //solverPush->writeDataToMeshObject(*mesh);
        //mesh->writeVTKFile( "out/solverPush.vtk" );
        //
        //solverSOA->writeDataToMeshObject(*mesh);
        //mesh->writeVTKFile( "out/solverSOA.vtk" );
        //
        //solverAOS->writeDataToMeshObject(*mesh);
        //mesh->writeVTKFile( "out/solverAOS.vtk" );

        // ====================================================================
        //              Output several files
        // ====================================================================


        //mesh->writeTimeSteps("out/timeSteps.dat");
        //mesh->writeTime("out/time.dat");
        //mesh->writeResultFields("out/ResultFields.dat");
        //mesh->writeOverviewFile("out/OverviewFile.dat");
        //mesh->writeConvergenceHistory("out/ConvergenceHistory.dat");

        // ========== Poiseuille Convergence Study ============================
        //ostringstream filename;
        //filename << "out/" << ny;
        //mesh->writeVelocityProfile(            ( filename.str() + "/VelocityProfile.dat" )    , 0.25);
        //mesh->writeResultFields(               ( filename.str() + "/ResultFields.dat" )            );
        //mesh->writeConvergenceHistory(         ( filename.str() + "/ConvergenceHistory.dat" )      );
        //mesh->writeOverviewFile(               ( filename.str() + "/OverviewFile.dat" )            );
        //mesh->writeVTKFile(                    ( filename.str() + "/ResultFields.vtk" )            );
        // ====================================================================

        // ========== Thermal Couette Convergence Study =======================
        //ostringstream filename;
        //filename << "out/ConvergenceStudyThermalCouette/" << ny;
        //mesh->writeTemperatureProfile(         ( filename.str() + "/TemperatureProfile.dat" )       , 0.5);
        //mesh->writeOverviewFile(               ( filename.str() + "/OverviewFile.dat" )                  );
        //mesh->writeVTKFile(                    ( filename.str() + "/ResultFields.vtk" )                  );
        // ====================================================================
        
        // ====================================================================
        //ostringstream filename;
        //filename << "out/Ra1e" << i + 3;
        //mesh->writeOverviewFile(  ( filename.str() + "/OverviewFile.dat" )     );
        //mesh->writeTimeSteps(     ( filename.str() + "/TimeSteps.dat"    )     );
        //mesh->writeVTKFile(       ( filename.str() + "/ResultFields.vtk" )     );
        //mesh->writeVelocityU(     ( filename.str() + "/VelocityU.dat" )        );
        //mesh->writeVelocityV(     ( filename.str() + "/VelocityV.dat" )        );
        // ====================================================================
        
        // ====================================================================
        //mesh->writeResultFields("out/ResultFields.dat");
        //mesh->writeOverviewFile("out/OverviewFile.dat");
        //mesh->writeConvergenceHistory("out/ConvergenceHistory.dat");
        ////mesh->writePressureGradientProfile("out/PressureGradientProfile.dat", 0.5);
        ////mesh->writeVelocityProfile("out/VelocityProfile.dat", 0.5);
        //mesh->writeTemperatureProfile("out/TemperatureProfile.dat", 0.5);
        //mesh->writeTimeSteps("out/TimeSteps.dat");
        ////mesh->writeVelocityU("out/VelocityU.dat");
        ////mesh->writeVelocityV("out/VelocityV.dat");
        ////mesh->writeTemperature("out/Temperature.dat");
        ////mesh->writeDensity("out/Density.dat");
        // ====================================================================

        system("pause");
        //delete solverPull;
        //delete solverPush;
        //delete solverSOA;
        //delete solverAOS;
    }
}