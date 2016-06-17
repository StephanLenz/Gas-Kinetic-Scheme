// ============================================================================
//
//                      Compressible Thermal GKS
//
// ============================================================================


#include "GKSMesh.h"
#include "BoundaryCondition.h"
#include <iostream>
#include <sstream>
#include <chrono>

using namespace std;

int main(int argc, char* argv[])
{
    /*

    // ========================================================================
    //
    //                  Couette-Flow (Periodic)
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 0.1;

    param.numberOfIterations = 1000000;
    param.outputIntervalVTK = 1000000;
    param.outputInterval = 10000;

    param.convergenceCriterium = 1.0e-8;

    param.CFL = 0.5;

    param.verbose = false;

    // ========================================================================

    FluidParameter fluidParam;

    fluidParam.K = 1;
    fluidParam.nu = 1e-2;
    fluidParam.R = 200.0;
    fluidParam.Force.x = 1e-4;
    fluidParam.Force.y = 0.0;

    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param, fluidParam);

    // Define Boundary Conditions
    //    -----------
    //    |    3    |
    //    | 0     2 |
    //    |    1    |
    //    -----------
    mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);

    // Generate Mesh
    mesh->generateRectMesh(W, H, 1, 64);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    */

    /*

    // ========================================================================
    //
    //                  Poiseuille-Flow (Force driven)
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 0.1;

    param.numberOfIterations = 100000;
    param.outputIntervalVTK = 1000;
    param.outputInterval = 1000;
    param.CFL = 0.5;
    
    param.convergenceCriterium = 1.0e-8;

    param.verbose = false;
    param.fluxOutput = true;

    // ========================================================================

    FluidParameter fluidParam;

    fluidParam.K  = 1;
    fluidParam.nu = 1e-2;
    fluidParam.R = 200.0;
    fluidParam.Force.x = 1e-4;
    fluidParam.Force.y = 0.0;

    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param, fluidParam);

    // Define Boundary Conditions
    //    -----------
    //    |    1    |
    //    |         |
    //    |    0    |
    //    -----------
    mesh->addInterfaceBoundaryCondition(0.0);
    mesh->addInterfaceBoundaryCondition(0.0);

    // Generate Mesh
    mesh->generateRectMeshPeriodicInterfaceBCs(W, H, 1, 64);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    */

    ///*

    // ========================================================================
    //
    //                  Driven Cavity
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 1000000;
    param.outputIntervalVTK = 1000000;
    param.outputInterval = 10000;

    param.convergenceCriterium = 1.0e-5;

    param.CFL = 0.1;
    param.fluxOutput = false;

    param.verbose = false;

    // ========================================================================

    FluidParameter fluidParam;

    fluidParam.K = 0;
    fluidParam.R = 208.0;

    double Re = 400;
    double Ma = 0.15;
    double T = 313.15;

    double cs = sqrt(fluidParam.R * T);
    double uTop = cs*Ma;

    fluidParam.nu = uTop / Re;

    double lambda = 1.0 / ( 2.0 * fluidParam.R * T );

    cout << "Re = " << uTop / fluidParam.nu << endl;
    cout << "Ma = " << uTop / sqrt(fluidParam.R * T) << endl;

    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param, fluidParam);

    // Define Boundary Conditions
    //    -----------
    //    |    3    |
    //    | 0     2 |
    //    |    1    |
    //    -----------
    mesh->addInterfaceBoundaryCondition(0.0);
    mesh->addInterfaceBoundaryCondition(0.0);
    mesh->addInterfaceBoundaryCondition(0.0);
    mesh->addInterfaceBoundaryCondition(uTop);

    // Generate Mesh
    mesh->generateRectMeshInterfaceBCs(W, H, 64, 64);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);

    //*/

    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================

    //cout << mesh->toString();

    mesh->writeMeshAsText("out/Mesh.txt");

    mesh->iterate();

    //mesh->writeTimeSteps("out/timeSteps.dat");

    mesh->writeVelocityU("out/VelocityU.dat");
    mesh->writeVelocityV("out/VelocityV.dat");

    //ostringstream filename;
    //filename << "out/PoiseuillePresGradConvergenceStudy/" << ny;
    //mesh->writeVelocityProfile(    ( filename.str() + "/VelocityProfile.dat" )  , 0.5);
    //mesh->writeConvergenceHistory( ( filename.str() + "/ConvergenceHistory.dat" )    );
    //mesh->writeOverviewFile(       ( filename.str() + "/OverviewFile.dat" )          );
    
    //ostringstream filename;
    //filename << "out/DrivenCavity/Re" << Re << "/" << nx;
    ////mesh->writeVelocityProfile(( filename.str() + "/VelocityProfile.dat" ), 0.5);
    //mesh->writeVelocityU(          (filename.str() + "/VelocityU.dat"          ));
    //mesh->writeVelocityV(          (filename.str() + "/VelocityV.dat"          ));
    //mesh->writeConvergenceHistory(( filename.str() + "/ConvergenceHistory.dat" ));
    //mesh->writeOverviewFile(      ( filename.str() + "/OverviewFile.dat" ));


    //mesh->writeConvergenceHistory("out/ConvergenceHistory.dat");
    mesh->writeOverviewFile("out/OverviewFile.dat");

    //char a; cin >> a;
}