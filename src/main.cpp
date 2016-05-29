// ============================================================================
//
//                      Compressible Thermal GKS
//
// ============================================================================


#include "GKSMesh.h"
#include "BoundaryCondition.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    /*

    // ========================================================================
    //
    //                  Couette-Flow
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 0.2;

    param.numberOfIterations = 1000000;
    param.outputInterval = 10000;
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
    mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.00, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.005, 0.0, 0.0);

    // Generate Mesh
    mesh->generateRectMeshPeriodic(W, H, 1, 128);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    */

    ///*

    // ========================================================================
    //
    //                  Poiseulle-Flow (Force driven)
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 0.2;

    param.numberOfIterations = 1000000;
    param.outputInterval = 10000;
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
    mesh->generateRectMeshPeriodic(W, H, 1, 128);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    //*/

    /*

    // ========================================================================
    //
    //                  Driven-Cavity
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 1000000;
    param.outputInterval = 100000;
    param.CFL = 0.5;

    param.verbose = false;

    // ========================================================================

    FluidParameter fluidParam;

    fluidParam.K = 1;
    fluidParam.nu = 1e-4;
    fluidParam.R = 200.0;

    double uTop = 0.01;
    double TAve = 10.0;

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
    mesh->generateRectMesh(W, H, 100, 100);

    // Initialize Values
    mesh->initMeshConstant(10.0, 0.0, 0.0, TAve);

    */

    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================

    mesh->applyBoundaryCondition();

    //cout << mesh->toString();

    mesh->iterate();

    //mesh->writeTimeSteps("out/timeSteps.dat");
    mesh->writeVelocityProfile("out/VelocityProfile.dat");
    
    //char a; cin >> a;
}