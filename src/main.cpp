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
    //double ReList[] = {40.0, 100.0, 400.0, 1000.0};
    //int    nyList[] = {8, 16, 32, 64, 128};
    //for (int i = 0; i < 4;i++)
    //for (int j = 0; j < 5;j++)
    {

        /*

        // ========================================================================
        //
        //                  Couette-Flow (Periodic)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 1000;
        param.outputIntervalVTK = 10;
        param.outputInterval = 10;

        param.convergenceCriterium = 1.0e-7;

        param.L = 1.0;
        param.CFL = 0.5;

        param.fluxOutput = false;
        param.resOutput = false;

        param.verbose = false;

        // ========================================================================

        FluidParameter fluidParam;

        fluidParam.K  = 1;
        fluidParam.nu = 1e-2;
        fluidParam.R = 208.0;
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
        mesh->addBoundaryCondition(1, 0, 0, 1,  1.0, 0.0,  0.0, 0.0);
        mesh->addBoundaryCondition(1, 0, 0, 1,  1.0, 0.1, 0.0, 0.0);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(incompressible, W, H, 1, 16);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, 3.0/2.0);

        */
    
        /*

        // ========================================================================
        //
        //                  Poiseuille-Flow (non periodic, pressure driven)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 100000;
        param.outputIntervalVTK = 10000;
        param.outputInterval = 10000;

        param.convergenceCriterium = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.1;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    ny = 32;
        double Re = 40;
        double u0 = 0.1;

        fluidParam.K = 1;
        fluidParam.nu = (u0*param.L)/Re;
        fluidParam.R = 208.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = 0.0;

        double dRho = 3.0 * (u0*8.0*fluidParam.nu) / (param.L*param.L);
    
        // ========== Diffusive Scaling ==========
        //int nyRef = 8;
        //int ny    = 256;
        //double Re = 40.0;

        //double uRef = 0.1;

        //double u    = uRef * ((double)nyRef/(double)ny);

        //fluidParam.K  = 1;
        //fluidParam.nu = u*param.L / Re;
        //fluidParam.R = 208.0;
        //fluidParam.Force.x = 8.0*u*fluidParam.nu / (param.L*param.L);
        //fluidParam.Force.y = 0.0;

        // ========== Test Paramesters ==========
        //fluidParam.K = 1;
        //fluidParam.nu = 1.0e-2;
        //fluidParam.R = 208.0;
        //fluidParam.Force.x = 1.0e-4;
        //fluidParam.Force.y = 0.0;

        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    3    |
        //    | 0     2 |
        //    |    1    |
        //    -----------
        mesh->addBoundaryCondition(0, 1, 1, 1, 1.0+dRho, 0.0, 0.0, 0.0);
        mesh->addBoundaryCondition(1, 0, 0, 1, 0.0     , 0.0, 0.0, 0.0);
        mesh->addBoundaryCondition(0, 1, 1, 1, 1.0     , 0.0, 0.0, 0.0);
        mesh->addBoundaryCondition(1, 0, 0, 1, 0.0     , 0.0, 0.0, 0.0);

        // Generate Mesh
        mesh->generateRectMesh(incompressible, W, H, ny, ny);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, 3.0/2.0);

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

        param.numberOfIterations = 10000000;
        param.outputIntervalVTK = 100000;
        param.outputInterval = 100000;

        param.convergenceCriterium = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.1;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    ny = 16;
        double Re = 40.0;
        double u0 = 0.1;

        fluidParam.K = 1;
        fluidParam.nu = (u0*param.L)/Re;
        fluidParam.R = 208.0;
        fluidParam.Force.x = (u0*8.0*fluidParam.nu) / (param.L*param.L);
        fluidParam.Force.y = 0.0;

        double T      = 293.15;
        double lambda = 1.0 / (2.0 * fluidParam.R * T);
    
        // ========== Diffusive Scaling ==========
        //int nyRef = 8;
        //int ny    = 256;
        //double Re = 40.0;

        //double uRef = 0.1;

        //double u    = uRef * ((double)nyRef/(double)ny);

        //fluidParam.K  = 1;
        //fluidParam.nu = u*param.L / Re;
        //fluidParam.R = 208.0;
        //fluidParam.Force.x = 8.0*u*fluidParam.nu / (param.L*param.L);
        //fluidParam.Force.y = 0.0;

        // ========== Test Paramesters ==========
        //fluidParam.K = 1;
        //fluidParam.nu = 1.0e-2;
        //fluidParam.R = 208.0;
        //fluidParam.Force.x = 1.0e-4;
        //fluidParam.Force.y = 0.0;

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
        mesh->generateRectMeshPeriodic(compressible, W, H, 1, ny);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);

        //*/

        /*

        // ========================================================================
        //
        //                  Poiseuille-Flow (Force driven, InterfaceBCs)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 10000000;
        param.outputIntervalVTK = 10000000;
        param.outputInterval = 100000;

        param.convergenceCriterium = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.3;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    ny = nyList[j];
        double Re = ReList[i];
        double u0 = 0.1;

        fluidParam.K = 1;
        fluidParam.nu = (u0*param.L)/Re;
        fluidParam.R = 208.0;
        fluidParam.Force.x = (u0*8.0*fluidParam.nu) / (param.L*param.L);
        fluidParam.Force.y = 0.0;
    
        // ========== Diffusive Scaling ==========
        //int nyRef = 8;
        //int ny    = 256;
        //double Re = 40.0;

        //double uRef = 0.1;

        //double u    = uRef * ((double)nyRef/(double)ny);

        //fluidParam.K  = 1;
        //fluidParam.nu = u*param.L / Re;
        //fluidParam.R = 208.0;
        //fluidParam.Force.x = 8.0*u*fluidParam.nu / (param.L*param.L);
        //fluidParam.Force.y = 0.0;

        // ========== Test Paramesters ==========
        //fluidParam.K = 1;
        //fluidParam.nu = 1.0e-2;
        //fluidParam.R = 208.0;
        //fluidParam.Force.x = 1.0e-4;
        //fluidParam.Force.y = 0.0;

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
        mesh->generateRectMeshPeriodicInterfaceBCs(incompressible, W, H, 1, ny);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, 3.0/2.0);

        */

        // ================================================================================================================================================
        // ================================================================================================================================================
        // ================================================================================================================================================
        // ================================================================================================================================================
        // ================================================================================================================================================

        //cout << mesh->toString();

        //mesh->writeMeshAsText("out/Mesh.txt");

        mesh->iterate();

        //mesh->writeTimeSteps("out/timeSteps.dat");

        //mesh->writeVelocityU("out/VelocityU.dat");
        //mesh->writeVelocityV("out/VelocityV.dat");

        //ostringstream filename;
        //filename << "out/PoiseuilleFlowIncompressible/ConvergenceStudyInterface/" << "Re" << Re << "/" << ny;
        //mesh->writeVelocityProfile(    ( filename.str() + "/VelocityProfile.dat" )  , 0.5);
        //mesh->writeConvergenceHistory( ( filename.str() + "/ConvergenceHistory.dat" )    );
        //mesh->writeOverviewFile(       ( filename.str() + "/OverviewFile.dat" )          );
        //mesh->writeVTKFile(            ( filename.str() + "/ResultFields.vtk" )          );
    
        //ostringstream filename;
        //filename << "out/DrivenCavity/Re" << Re << "/" << nx;
        ////mesh->writeVelocityProfile(( filename.str() + "/VelocityProfile.dat" ), 0.5);
        //mesh->writeVelocityU(          (filename.str() + "/VelocityU.dat"          ));
        //mesh->writeVelocityV(          (filename.str() + "/VelocityV.dat"          ));
        //mesh->writeConvergenceHistory(( filename.str() + "/ConvergenceHistory.dat" ));
        //mesh->writeOverviewFile(      ( filename.str() + "/OverviewFile.dat" ));


        mesh->writeConvergenceHistory("out/ConvergenceHistory.dat");
        mesh->writeOverviewFile("out/OverviewFile.dat");

        //char a; cin >> a;
    }
}