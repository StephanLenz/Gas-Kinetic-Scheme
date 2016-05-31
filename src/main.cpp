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
    ///*

    // ========================================================================
    //
    //                  Couette-Flow (Periodic)
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 100000;
    param.outputInterval = 100000;
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
    mesh->generateRectMeshPeriodicInterfaceBCs(W, H, 3, 3);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    //*/
    
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================

    //cout << mesh->toString();

    mesh->writeMeshAsText("out/Mesh.txt");

    //esh->iterate();

    //mesh->writeTimeSteps("out/timeSteps.dat");
    //mesh->writeVelocityProfile("out/VelocityProfile.dat");
    
    //char a; cin >> a;
}