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
    double W = 0.25;

    param.numberOfIterations = 10000;
    param.outputInterval = 10000;
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
    double W = 0.1;

    param.numberOfIterations = 100000;
    param.outputInterval = 100;
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
    int ny = 64;
    mesh->generateRectMeshPeriodicInterfaceBCs(W, H, 6, ny);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, 1.0);

    //*/
    
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================

    //cout << mesh->toString();

    //mesh->writeMeshAsText("out/Mesh.txt");

    mesh->iterate();

    ostringstream filename;
    filename << "out/VelocityProfileForcingInterfaceBC" << ny << ".dat";
    mesh->writeVelocityProfile(filename.str(), 0.05);
    
    //char a; cin >> a;
}