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

#include "GKSSolver.h"
#include "GKSMesh.h"
#include "BoundaryCondition.h"
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
    {
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
        param.ghostOutput = false;                      // include ghost cells in VTK files
        param.csvOutput   = true;                       // output csv files for postprocessing

        param.numberOfIterations = 10000;         // maximal number of Iterations
        param.outputIntervalVTK = 10000;                 // Output interval for VTK Files (and .dat files)
        param.outputInterval = 1000;                    // Output interval for Output on the screen

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

        // domain size in [m] for mesh generation
        double H = 1.0;
        double W = 1.0;

        FluidParameter fluidParam;

        int    nx = 32;          // number of cells in x direction
        int    ny = 32;         // number of cells in y direction
        double Re = 4.0;        // Reynolds number
        double u0 = 0.1;        // Velocits in the mid of the channel
        param.L = 1.0;          // reference length for Re number

        fluidParam.K = 1.0;                                                   // internal degrees of freedom
        fluidParam.nu = (u0*param.L)/Re;                                    // kinematic viskosity
        fluidParam.R = 200.0;                                               // specific gas constant
        fluidParam.Force.x = (u0*8.0*fluidParam.nu) / (param.L*param.L);    // acceleration in x direction [m/s^2]
        fluidParam.Force.y = 0.0;                                           // acceleration in y direction [m/s^2]
        fluidParam.BoussinesqForce.x = 0.0;                                 // acceleration only allpied to density variations [m/s^2]
        fluidParam.BoussinesqForce.y = 0.0;                                 // acceleration only allpied to density variations [m/s^2]
        fluidParam.rhoReference = 1.0;                                      // reference density
        fluidParam.Pr = 1.0;                                                // Prandl number 

        double T      = 1.0;                                              // reference temperature [K]
        double lambda = 1.0 / (2.0 * fluidParam.R * T);
        
        // ========================================================================
        //                  Definition of boundary conditions
        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    3    |
        //    | 0     2 |
        //    |    1    |
        //    -----------
        mesh->addBoundaryCondition(periodic, 0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(wall, 0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(periodic, 0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(wall, 0.0, 0.0, 0.0, lambda);

        // Generate the Mesh
        mesh->generateRectMeshGraded(compressible, W, H, nx, ny, 1.0, 1.0);

        // Set initial condition
        mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);
        //mesh->initMeshParabularVelocity(1.0, u0, 0.0, lambda);
        //mesh->initMeshSineVelocity(1.0, u0, 0.0, lambda);

        //*/
        
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
        //mesh->iterate();
        // ====================================================================
        // ====================================================================

        GKSSolver* solver = new GKSSolver(param, fluidParam);

        solver->readMeshFromMeshObject(*mesh);

        mesh->iterate();

        solver->iterate();

        //solver->writeDataToMeshObject(*mesh);

        //mesh->writeVTKFile( "out/solver.vtk" );

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

        //system("pause");
        delete mesh;
        delete solver;
    }
}