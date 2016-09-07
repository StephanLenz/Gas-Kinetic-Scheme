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

# define M_PI           3.14159265358979323846  /* pi */

using namespace std;

int main(int argc, char* argv[])
{
    double ReList[] = {40.0, 100.0, 400.0, 1000.0};
    int    nyList[] = {8, 16, 32, 64};
    double RaList[] = {1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8};
    //for (int i = 0; i < 4;i++)      //ReList
    //for (int j = 0; j < 4;j++)      // nyList
    //for(int i = 0; i < 6; i++)      // RaList
    {
    
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
        param.outputIntervalVTK = 10000;
        param.outputInterval = 10000;

        param.convergenceCriterium[0] = 1.0;
        param.convergenceCriterium[1] = -1.0e-10;
        param.convergenceCriterium[2] = 1.0;
        param.convergenceCriterium[3] = 1.0;

        param.CFL = 0.5;

        param.verbose = false;
        param.fluxOutput = true;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    nx = 2;
        int    ny = 2;//nyList[j];
        double Re = 40.0;
        double u0 = 0.1;
        double angle = atan(0.0);
        param.L = 1.0*cos(angle);

        fluidParam.K = 1;
        fluidParam.nu = (u0*param.L)/Re;
        fluidParam.R = 200.0;
        fluidParam.Force.x = 0.0;//cos(angle) * (u0*8.0*fluidParam.nu) / (param.L*param.L);
        fluidParam.Force.y = 0.0;//sin(angle) * (u0*8.0*fluidParam.nu) / (param.L*param.L);
        fluidParam.BoussinesqForce.x = 0.0;
        fluidParam.BoussinesqForce.y = 0.0;
        fluidParam.rhoReference = 1.0;

        double T      = 300.0;
        double lambda = 1.0 / (2.0 * fluidParam.R * T);

        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    3    |
        //    | 0     2 |
        //    |    1    |
        //    -----------
        mesh->addBoundaryCondition(periodic, 0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(periodic, 0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(periodic, 0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(periodic, 0.0, 0.0, 0.0, lambda);

        Interface::setInterpolationOrder(1);

        // Generate Mesh
        mesh->generateRectMeshGraded(compressible, W, H, nx, ny, 1.0, 1.0);

        // Initialize Values
        //mesh->initMeshConstant(1.0, 1.0, 0.0, lambda);
        //mesh->initMeshParabularVelocity(1.0, u0, 0.0, lambda);
        mesh->initMeshSineVelocity(1.0, 1.0, 0.0, lambda);

        //double rhoLin[] = {1.0, 1.0};
        //double uLin[] = {0.0, 1.0};
        //double vLin[] = {0.0, 0.0};
        //double lambdaLin[] = {lambda, lambda};
        //mesh->initMeshLinear(rhoLin, uLin, vLin, lambdaLin);

        */


        /*

        // ========================================================================
        //
        //                  Thermal Couette-Flow (Periodic)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 1000000;
        param.outputIntervalVTK = 100000;
        param.outputInterval = 100000;

        param.convergenceCriterium[0] = 1.0e-10;
        param.convergenceCriterium[1] = 1.0e-10;
        param.convergenceCriterium[2] = 1.0e-8;
        param.convergenceCriterium[3] = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.5;

        param.fluxOutput = true;
        param.resOutput = false;

        param.verbose = false;

        // ========================================================================

        FluidParameter fluidParam;

        int    ny = nyList[j];
        int    nx = 1;
        double Re = 1.0;//ReList[i];
        double u0 = 1.0;

        fluidParam.K = 1;
        fluidParam.nu = 0.1;
        fluidParam.R = 200.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = 0.0;

        double TTop   = 293.15 + 5.0e-4;
        double TBot   = 293.15;
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
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[0]);
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, u0 , 0.0, lambda[1]);
        //mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0);
        //mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, u0 , 0.0, 0.0);

        Interface::setInterpolationOrder(1);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(compressible, W, H, 1, ny);

        // Initialize Values
        //mesh->initMeshConstant(1.0, 0.0, 0.0, 1.5 );
        mesh->initMeshLinear(rho, U, V, lambda);

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
        param.outputIntervalVTK = 1000;
        param.outputInterval = 1000;

        param.convergenceCriterium[0] = 1.0;
        param.convergenceCriterium[1] = 1.0e-10;
        param.convergenceCriterium[2] = 1.0;
        param.convergenceCriterium[3] = 1.0;

        param.CFL = 0.7;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    nx = 2;
        int    ny = 2;//nyList[j];
        double Re = 4.0;
        double u0 = 0.1;
        double angle = 0.0*M_PI;//atan(0.0);
        param.L = 1.0;//*cos(angle);

        fluidParam.K = 1;
        fluidParam.nu = (u0*param.L)/Re;
        fluidParam.R = 200.0;
        fluidParam.Force.x = cos(angle) * (u0*8.0*fluidParam.nu) / (param.L*param.L);
        fluidParam.Force.y = sin(angle) * (u0*8.0*fluidParam.nu) / (param.L*param.L);
        fluidParam.BoussinesqForce.x = 0.0;
        fluidParam.BoussinesqForce.y = 0.0;
        fluidParam.rhoReference = 1.0;

        double T      = 300.0;
        //double lambda = 1.0 / (2.0 * fluidParam.R * T);
        double lambda = 1.5;

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
        //mesh->addBoundaryCondition(periodic, 0.0, 0.0, 0.0, lambda);

        Interface::setInterpolationOrder(1);

        // Generate Mesh
        mesh->generateRectMeshGraded(incompressible, W, H, nx, ny, 1.0, 1.0);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);
        //mesh->initMeshParabularVelocity(1.0, u0, 0.0, lambda);
        //mesh->initMeshSineVelocity(1.0, u0, 0.0, lambda);

        //*/
    
        /*

        // ========================================================================
        //
        //                  Atmospheric Pressure
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 10000000;
        param.outputIntervalVTK = 100000;
        param.outputInterval = 100000;

        param.convergenceCriterium[0] = -1.0e-10;
        param.convergenceCriterium[1] = -1.0e-10;
        param.convergenceCriterium[2] = -1.0e-10;
        param.convergenceCriterium[3] = -1.0e-10;

        param.L = 1.0;
        param.CFL = 0.7;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        int    nx = 1;
        int    ny = 32;//nyList[j];

        double Ra = 1000;

        double TReference = 300.0;
        double TTop   = TReference - 0.0;
        double TBot   = TReference + 0.0;
        double g      = 10.0;

        fluidParam.K = 1;
        fluidParam.nu = sqrt( (g * H)/Ra * (TBot - TTop)/TReference );
        fluidParam.R = 200.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = -g;
        fluidParam.BoussinesqForce.x = 0.0;
        fluidParam.BoussinesqForce.y = 0.0;
        fluidParam.rhoReference = 1.0;

        double lambdaReference = 1.0 / (2.0 * fluidParam.R * TReference);
        double lambda[] = { 1.0 / (2.0 * fluidParam.R * TBot), 1.0 / (2.0 * fluidParam.R * TTop) };
        double rho[]    = { fluidParam.rhoReference * lambda[0]/lambdaReference, fluidParam.rhoReference * lambda[1]/lambdaReference };
        double U[] = { 0.0, 0.0 };
        double V[] = { 0.0, 0.0 };

        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    3    |
        //    | 0     2 |
        //    |    1    |
        //    -----------
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0,    0.0, 0.0, lambda[0]);
        //mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0   );
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0,    0.0, 0.0, lambda[1]);
        //mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0   );

        Interface::setInterpolationOrder(1);
        
        // Generate Mesh
        //mesh->generateRectMesh(compressible, W, H, nx, ny);
        mesh->generateRectMeshPeriodic(compressible, W, H, nx, ny);

         // Initialize Values
        //mesh->initMeshConstant(1.0, 0.0, 0.0, lambda[0]);
        //mesh->initMeshLinear(rho, U, V, lambda);
        mesh->initMeshAtmospheric(fluidParam.rhoReference, 0.0, 0.0, lambdaReference, g);

        */
    
        /*

        // ========================================================================
        //
        //                  Thermal Driven Cavity
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 10000;
        param.outputIntervalVTK = 1000;
        param.outputInterval = 1000;

        param.convergenceCriterium[0] = 1.0e-10;
        param.convergenceCriterium[1] = 1.0e-10;
        param.convergenceCriterium[2] = 1.0e-10;
        param.convergenceCriterium[3] = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.3;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        int    nx = 16;
        int    ny = 16;//nyList[j];

        double Ra = 1.0e5;

        double TReference = 300.0;
        double TTop   = TReference - 100.0;
        double TBot   = TReference + 100.0;
        double g      = 10.0;

        fluidParam.K = 1;
        fluidParam.nu = sqrt( (g * H)/Ra * (TBot - TTop)/TReference );
        fluidParam.R = 200.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = -g;
        fluidParam.BoussinesqForce.x = 0.0;
        fluidParam.BoussinesqForce.y = 0.0;
        fluidParam.rhoReference = 1.0;

        double lambdaReference = 1.0 / (2.0 * fluidParam.R * TReference);
        double lambda[] = { 1.0 / (2.0 * fluidParam.R * TBot), 1.0 / (2.0 * fluidParam.R * TTop) };
        double rho[]    = { fluidParam.rhoReference * lambda[0]/lambdaReference, fluidParam.rhoReference * lambda[1]/lambdaReference };
        double U[] = { 0.0, 0.0 };
        double V[] = { 0.0, 0.0 };

        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    3    |
        //    | 0     2 |
        //    |    1    |
        //    -----------
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[0]);
        mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0   );
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[1]);
        mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0   );

        Interface::setInterpolationOrder(1);
        
        // Generate Mesh
        //mesh->generateRectMesh(compressible, W, H, nx, ny);
        //mesh->generateRectMesh(compressible, W, H, nx, ny);
        mesh->generateRectMeshGraded(compressible, W, H, nx, ny, 0.1, 0.1);

         // Initialize Values
        //mesh->initMeshConstant(1.0, 0.0, 0.0, 0.5*(lambda[0] + lambda[1]));
        mesh->initMeshLinearHorizontal(rho, U, V, lambda);
        //mesh->initMeshParabularVelocity(1.0, u0, 0.0, lambda);

        */
    
        /*

        // ========================================================================
        //
        //                  Rayleigh-Bernard Convection
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 10000000;
        param.outputIntervalVTK = 10000;
        param.outputInterval = 10000;

        param.convergenceCriterium[0] = 1.0e-10;
        param.convergenceCriterium[1] = 1.0e-10;
        param.convergenceCriterium[2] = 1.0e-10;
        param.convergenceCriterium[3] = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.5;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        int    nx = 32;
        int    ny = 32;//nyList[j];

        double Ra = 10000.0;

        double TReference = 300.0;
        double TTop   = TReference - 100.0;
        double TBot   = TReference + 100.0;
        double g      = 10.0;

        fluidParam.K = 1;
        fluidParam.nu = sqrt( (g * H*H*H)/Ra * (TBot - TTop)/TReference );
        fluidParam.R = 200.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = -g;
        fluidParam.BoussinesqForce.x = 0.0;
        fluidParam.BoussinesqForce.y = 0.0;
        fluidParam.rhoReference = 1.0;

        double lambdaReference = 1.0 / (2.0 * fluidParam.R * TReference);
        double lambda[] = { 1.0 / (2.0 * fluidParam.R * TBot), 1.0 / (2.0 * fluidParam.R * TTop) };
        double rho[]    = { fluidParam.rhoReference * lambda[0]/lambdaReference, fluidParam.rhoReference * lambda[1]/lambdaReference };
        double U[] = { 0.0, 0.0 };
        double V[] = { 0.0, 0.0 };

        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------
        //    |    3    |
        //    | 0     2 |
        //    |    1    |
        //    -----------
        mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0   );
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[0]);
        mesh->addBoundaryCondition(1, 0, 0, 1,  0.0, 0.0, 0.0, 0.0   );
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[1]);

        Interface::setInterpolationOrder(1);
        
        // Generate Mesh
        //mesh->generateRectMesh(compressible, W, H, nx, ny);
        mesh->generateRectMesh(compressible, W, H, nx, ny);

         // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, 0.5*(lambda[0] + lambda[1]));
        //mesh->initMeshLinear(rho, U, V, lambda);
        //mesh->initMeshParabularVelocity(1.0, u0, 0.0, lambda);

        */

        // ================================================================================================================================================
        // ================================================================================================================================================
        // ================================================================================================================================================
        // ================================================================================================================================================
        // ================================================================================================================================================

        //cout << mesh->toString();

        //mesh->writeMeshAsText("out/Mesh.txt");

        mesh->writeVTKFile("out/InitialState.vtk");
        mesh->writeVTKFileFlux("out/InitialStateFlux.vtk");

        mesh->iterate();

        //mesh->writeTimeSteps("out/timeSteps.dat");

        // ========== Poiseuille Convergence Study ============================
        //ostringstream filename;
        //filename << "out/" << ny;
        //mesh->writeVelocityProfile(            ( filename.str() + "/VelocityProfile.dat" )    , 0.5);
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
        mesh->writeResultFields("out/ResultFields.dat");
        mesh->writeOverviewFile("out/OverviewFile.dat");
        //mesh->writeConvergenceHistory("out/ConvergenceHistory.dat");
        ////mesh->writePressureGradientProfile("out/PressureGradientProfile.dat", 0.5);
        ////mesh->writeVelocityProfile("out/VelocityProfile.dat", 0.5);
        ////mesh->writeTemperatureProfile("out/TemperatureProfile.dat", 0.5);
        //mesh->writeTimeSteps("out/TimeSteps.dat");
        ////mesh->writeVelocityU("out/VelocityU.dat");
        ////mesh->writeVelocityV("out/VelocityV.dat");
        ////mesh->writeTemperature("out/Temperature.dat");
        ////mesh->writeDensity("out/Density.dat");
        // ====================================================================

        //system("pause");
    }
}