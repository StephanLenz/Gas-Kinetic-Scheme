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
    ///*

    // ========================================================================
    //
    //                  Couette-Flow (Periodic)
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 0.25;

    param.numberOfIterations = 1000000;
    param.outputIntervalVTK = 1000000;
    param.outputInterval = 10000;

    param.convergenceCriterium = 1.0e-6;

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
    mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.00, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.00, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.00, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.005, 0.0, 0.0);

    // Generate Mesh
    mesh->generateRectMesh(W, H, 2, 2);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    */

    ///*

    // ========================================================================
    //
    //                  Poiseuille-Flow (Force driven)
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 100;
    param.outputInterval = 1;
    param.CFL = 0.5;

    param.verbose = false;

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

    //*/
    
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param, fluidParam);

    // Define Boundary Conditions
    //    -----------
    //    |    3    |
    //    | 0     2 |
    //    |    1    |
    //    -----------
    mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.0, 0.0, 1.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 1.0);
    mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.0, 0.0, 1.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 1.0);

    // Generate Mesh
    mesh->generateRectMesh(W, H, 100, 100);

    // Initialize Values
    // mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    double rho[] = { 1.0, 1.0 + 1.0e-3 };

    mesh->initMeshLinearDensity(rho, 0.0, 0.0, lambda);

    */

    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================

    //cout << mesh->toString();

    mesh->writeMeshAsText("out/Mesh.txt");

    //esh->iterate();

    //mesh->writeTimeSteps("out/timeSteps.dat");

    //mesh->writeVelocityU("out/VelocityU.dat");
    //mesh->writeVelocityV("out/VelocityV.dat");

    //ostringstream filename;
    //filename << "out/PoiseuillePresGradConvergenceStudy/" << ny;
    //mesh->writeVelocityProfile(    ( filename.str() + "/VelocityProfile.dat" )  , 0.5);
    //mesh->writeConvergenceHistory( ( filename.str() + "/ConvergenceHistory.dat" )    );
    //mesh->writeOverviewFile(       ( filename.str() + "/OverviewFile.dat" )          );
    
    ostringstream filename;
    filename << "out/DrivenCavity/Re" << Re << "/" << nx;
    //mesh->writeVelocityProfile(( filename.str() + "/VelocityProfile.dat" ), 0.5);
    mesh->writeVelocityU(          (filename.str() + "/VelocityU.dat"          ));
    mesh->writeVelocityV(          (filename.str() + "/VelocityV.dat"          ));
    mesh->writeConvergenceHistory(( filename.str() + "/ConvergenceHistory.dat" ));
    mesh->writeOverviewFile(      ( filename.str() + "/OverviewFile.dat" ));


    //mesh->writeConvergenceHistory("out/ConvergenceHistory.dat");
    //mesh->writeOverviewFile("out/OverviewFile.dat");

    //char a; cin >> a;
}