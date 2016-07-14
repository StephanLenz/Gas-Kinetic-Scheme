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
    int    nyList[] = {8, 16, 32, 64};
    //for (int i = 0; i < 4;i++)
    //for (int j = 0; j < 4;j++)
    {


        /*

        // ========================================================================
        //
        //                  Thermal Couette-Flow (Periodic)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 10;
        param.outputIntervalVTK = 1;
        param.outputInterval = 1;

        param.convergenceCriterium[0] = 1.0e-10;
        param.convergenceCriterium[1] = 1.0e-10;
        param.convergenceCriterium[2] = 1.0e-10;
        param.convergenceCriterium[3] = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.5;

        param.fluxOutput = true;
        param.resOutput = false;

        param.verbose = false;

        // ========================================================================

        FluidParameter fluidParam;

        int    ny = 3;//nyList[j];
        int    nx = 1;
        double Re = 1.0;//ReList[i];
        double u0 = 1.0;

        fluidParam.K = 1;
        fluidParam.nu = 0.1;
        fluidParam.R = 200.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = 0.0;

        double TTop   = 10.0 + 5.0e-5;
        double TBot   = 10.0;
        double lambda[] = { 1.0 / (2.0 * fluidParam.R * TBot), 1.0 / (2.0 * fluidParam.R * TTop) };
        double rho[]    = { 1.0, 1.0 * lambda[1] / lambda[0] };
        //double rho[]    = { 1.0, 1.0 };
        double U[] = { 0.0, u0 };
        double V[] = { 0.0, 0.0 };


        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    1    |
        //    |         |
        //    |    0    |
        //    -----------
        //mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[0]);
        //mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, u0 , 0.0, lambda[1]);
        mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0);
        mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, u0 , 0.0, 0.0);

        Interface::setInterpolationOrder(1);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(compressible, W, H, 1, ny);

        // Initialize Values
        //mesh->initMeshConstant(1.0, 0.0, 0.0, 1.5 );
        //mesh->initMeshLinear(rho, U, V, lambda);

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

        param.convergenceCriterium[0] = 1.0;
        param.convergenceCriterium[1] = 1.0e-10;
        param.convergenceCriterium[2] = 1.0;
        param.convergenceCriterium[3] = 1.0;

        param.L = 1.0;
        param.CFL = 0.7;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    nx = 1;
        int    ny = 16;//nyList[j];
        double Re = 40.0;
        double u0 = 0.1;

        fluidParam.K = 1;
        fluidParam.nu = (u0*param.L)/Re;
        fluidParam.R = 208.0;
        fluidParam.Force.x = (u0*8.0*fluidParam.nu) / (param.L*param.L);
        fluidParam.Force.y = 0.0;

        double T      = 293.15;
        double lambda = 1.0 / (2.0 * fluidParam.R * T);

        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    1    |
        //    |         |
        //    |    0    |
        //    -----------
        mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0 , 0.0, lambda);

        Interface::setInterpolationOrder(1);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(compressible, W, H, nx, ny);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);
        //mesh->initMeshParabularVelocity(1.0, u0, 0.0, lambda);

        //*/
    
        /*

        // ========================================================================
        //
        //                  Uniform Advection
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 1000000;
        param.outputIntervalVTK = 1000;
        param.outputInterval = 1000;

        param.convergenceCriterium[0] = 1.0e-10;
        param.convergenceCriterium[1] = -1.0e-10;
        param.convergenceCriterium[2] = 1.0e-10;
        param.convergenceCriterium[3] = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.7;

        param.verbose = false;
        param.fluxOutput = true;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    nx = 1;
        int    ny = 16;
        double Re = 40.0;
        double u0 = 0.1;
        double v0 = 0.0;

        fluidParam.K = 1;
        fluidParam.nu = (u0*param.L)/Re;
        fluidParam.R = 208.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = 0.0;

        double T      = 293.15;
        double lambda = 1.0 / (2.0 * fluidParam.R * T);

        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    1    |
        //    |         |
        //    |    0    |
        //    -----------

        Interface::setInterpolationOrder(1);

        // Generate Mesh
        mesh->generateRectMeshPeriodicTwoDirections(compressible, W, H, nx, ny);

        // Initialize Values
        //mesh->initMeshConstant(1.0, u0, v0, 1.5);
        //mesh->initMeshParabularVelocity(1.0, u0, 0.0, lambda);
        mesh->initMeshSineVelocity(1.0, u0, 0.0, lambda);

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
        //filename << "out/PoiseuilleFlowCompressible/ConvergenceStudy1stOrder/" << ny;
        //mesh->writeVelocityProfile(            ( filename.str() + "/VelocityProfile.dat" )          , 0.5);
        //mesh->writePressureGradientProfile(    ( filename.str() + "/PressureGradientProfile.dat" )  , 0.5);
        //mesh->writeConvergenceHistory(         ( filename.str() + "/ConvergenceHistory.dat" )            );
        //mesh->writeOverviewFile(               ( filename.str() + "/OverviewFile.dat" )                  );
        //mesh->writeVTKFile(                    ( filename.str() + "/ResultFields.vtk" )                  );


        mesh->writeConvergenceHistory("out/ConvergenceHistory.dat");
        mesh->writeOverviewFile("out/OverviewFile.dat");
        //mesh->writePressureGradientProfile("out/PressureGradientProfile.dat", 0.5);
        //mesh->writeVelocityProfile("out/VelocityProfile.dat", 0.5);
        mesh->writeTimeSteps("out/TimeSteps.dat");

        //system("pause");
    }
}