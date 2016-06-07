// ============================================================================
//
//                      Compressible Thermal GKS
//
// ============================================================================


#include "GKSMesh.h"
#include "BoundaryCondition.h"
#include <iostream>
#include <sstream>

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
    double W = 0.5;

    param.numberOfIterations = 10000;
    param.outputInterval = 10;
    param.CFL = 0.5;

    param.verbose = false;

    // ========================================================================

    FluidParameter fluidParam;

    fluidParam.K  = 1;
    fluidParam.nu = 1e-2;
    fluidParam.R = 200.0;
    fluidParam.Force.x = 0.0;
    fluidParam.Force.y = 0.0;

    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param, fluidParam);

    // Define Boundary Conditions
    //    -----------
    //    |    1    |
    //    |         |
    //    |    0    |
    //    -----------
    mesh->addBoundaryCondition(1, 1, 1, 1,  1.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 1, 1, 1,  1.0, 0.0, 0.0, 0.0);

    // Generate Mesh
    mesh->generateRectMeshPeriodic(W, H, 1, 8);

    // Initialize Values
    //mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    double rho[] = { 1.0, 1.0 + 1.0e-3 };

    mesh->initMeshLinearDensity(rho, 0.0, 0.0, 1.0);

    */
    
    ///*

    // ========================================================================
    //
    //                  Poiseuille-Flow (non periodic)
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 500000;
    param.outputInterval = 500000;
    param.CFL = 0.5;

    param.verbose = false;

    // ========================================================================

    FluidParameter fluidParam;

    fluidParam.K = 1;
    fluidParam.nu = 1e-2;
    fluidParam.R = 200.0;
    fluidParam.Force.x = 0.0;
    fluidParam.Force.y = 0.0;

    double dp   = 1.0e-4;
    double lambda = 1.0;
    double drho = 2.0*dp*lambda;
    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param, fluidParam);

    // Define Boundary Conditions
    //    -----------
    //    |    3    |
    //    | 0     2 |
    //    |    1    |
    //    -----------
    //mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0);
    //mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0);
    //mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0);
    //mesh->addBoundaryCondition(1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0);

    mesh->addBoundaryCondition(0, 1, 1, 1, 1.0+drho, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0     , 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(0, 1, 1, 1, 1.0     , 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0     , 0.0, 0.0, 0.0);

    // Generate Mesh
    int ny = 512;
    mesh->generateRectMesh(W, H, 33, ny);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);

    //double rho[] = { 1.0, 1.0 + 1.0e-3 };

    //mesh->initMeshLinearDensity(rho, 0.0, 0.0, 1.0);

    //*/

    /*

    // ========================================================================
    //
    //                  Poiseuille-Flow (Force driven)
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 1000000;
    param.outputInterval = 1000000;
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
    mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0);

    // Generate Mesh
    int ny = 512;
    mesh->generateRectMeshPeriodic(W, H, 1, ny);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    */

    /*

    // ========================================================================
    //
    //                  Driven-Cavity
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 100000;
    param.outputInterval = 100000;
    param.CFL = 0.1;

    param.verbose = false;

    // ========================================================================

    FluidParameter fluidParam;

    fluidParam.K = 1;
    fluidParam.nu = 1e-2;
    fluidParam.R = 200.0;

    double uTop = 0.1;

    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param, fluidParam);

    // Define Boundary Conditions
    //    -----------
    //    |    3    |
    //    | 0     2 |
    //    |    1    |
    //    -----------
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, uTop, 0.0, 0.0);

    // Generate Mesh
    mesh->generateRectMesh(W, H, 32, 32);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);
    
    */

    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================

    //cout << mesh->toString();

    //mesh->writeMeshAsText("out/Mesh.txt");

    mesh->iterate();

    //mesh->writeTimeSteps("out/timeSteps.dat");

    ostringstream filename;
    filename << "out/VelocityProfilePresGrad" << ny << ".dat";
    mesh->writeVelocityProfile(filename.str(), 0.5);
    
    //char a; cin >> a;
}