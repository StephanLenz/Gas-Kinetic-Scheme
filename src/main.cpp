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
    int    nyList[] = {8, 16, 32, 64, 128};
    //for (int i = 0; i < 4;i++)
    //for (int j = 0; j < 5;j++)
    {


        ///*

        // ========================================================================
        //
        //                  Thermal Couette-Flow (Periodic)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 1000000;
        param.outputIntervalVTK = 10000;
        param.outputInterval = 10000;

        param.convergenceCriterium = 1.0e-7;

        param.L = 1.0;
        param.CFL = 0.1;

        param.fluxOutput = false;
        param.resOutput = false;

        param.verbose = false;

        // ========================================================================

        FluidParameter fluidParam;

        int    ny = 32;//nyList[j];
        int    nx = 1;
        double Re = 40.0;//ReList[i];
        double u0 = 1.0;

        fluidParam.K = 1;
        fluidParam.nu = (u0*param.L)/Re;
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
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[0]);
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, u0 , 0.0, lambda[1]);

        Interface::setInterpolationOrder(3);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(compressible, W, H, 1, ny);

        // Initialize Values
        //mesh->initMeshConstant(1.0, 0.0, 0.0, 0.5 * (lambda[0] + lambda[1]) );
        mesh->initMeshLinear(rho, U, V, lambda);

        //*/
    
        /*

        // ========================================================================
        //
        //                  Poiseuille-Flow (non periodic, pressure driven)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 2.0;

        param.numberOfIterations = 100000000;
        param.outputIntervalVTK = 10000000;
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

        double dRho = W * 3.0 * (u0*8.0*fluidParam.nu) / (param.L*param.L);
    
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
        mesh->addBoundaryCondition(0, 1, 1, 1,    1.0+dRho, 0.0, 0.0, 0.0);
        mesh->addBoundaryCondition(1, 0, 0, 1,    0.0 ,     0.0, 0.0, 0.0);
        mesh->addBoundaryCondition(0, 1, 1, 1,    1.0 ,     0.0, 0.0, 0.0);
        mesh->addBoundaryCondition(1, 0, 0, 1,    0.0 ,     0.0, 0.0, 0.0);

        Interface::setInterpolationOrder(3);

        // Generate Mesh
        mesh->generateRectMesh(incompressible, W, H, 2*ny+1, ny+1);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, 3.0/2.0);

        */

        /*

        // ========================================================================
        //
        //                  Poiseuille-Flow (Force driven)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 100000000;
        param.outputIntervalVTK = 1000000;
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
        int    ny = 32;//nyList[j];
        double Re = 40.0;//ReList[i];
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

        Interface::setInterpolationOrder(3);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(compressible, W, H, 1, ny);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);

        */

        /*

        // ========================================================================
        //
        //                  Poiseuille-Flow (Force driven)
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 100000000;
        param.outputIntervalVTK = 1000000;
        param.outputInterval = 10000;

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
        double Re = 40.0;//ReList[i];
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

        Interface::setInterpolationOrder(1);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(incompressible, W, H, 1, ny);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, 3.0/2.0);

        */

        /*

        // ========================================================================
        //
        //                  Pressure Wave
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 18000;
        param.outputIntervalVTK = 60;
        param.outputInterval = 60;

        param.convergenceCriterium = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.05;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    ny = 32;//nyList[j];
        double Re = 40.0;//ReList[i];
        double u0 = 0.1;

        fluidParam.K = 1;
        fluidParam.nu = 0.01;
        fluidParam.R = 208.0;
        fluidParam.Force.x = 0.0;
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

        Interface::setInterpolationOrder(3);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(compressible, W, H, ny, ny);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);
        double rhoValues[] = {10.0, 9.0};
        mesh->initMeshLinearDensity(rhoValues, 0.0, 0.0, lambda);

        */

        /*

        // ========================================================================
        //
        //                  Hydrostatic Pressure
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 100000;
        param.outputIntervalVTK = 100;
        param.outputInterval = 100;

        param.convergenceCriterium = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.05;

        param.verbose = false;
        param.fluxOutput = false;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Weidongs Parameters ==========
        int    ny = 32;//nyList[j];
        double Re = 40.0;//ReList[i];
        double u0 = 0.1;

        fluidParam.K = 1;
        fluidParam.nu = 0.1;
        fluidParam.R = 208.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = -9.81;

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
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda);
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda);

        Interface::setInterpolationOrder(3);

        // Generate Mesh
        mesh->generateRectMesh(compressible, W, H, ny, ny);

        // Initialize Values
        mesh->initMeshConstant(1.0, 0.0, 0.0, lambda);
        //mesh->initMeshLinearDensity(1.0, 0.0, 0.0, lambda);

        */

        /*

        // ========================================================================
        //
        //                  Rayleigh-Bernard
        //
        // ========================================================================

        Parameters param;

        double H = 1.0;
        double W = 1.0;

        param.numberOfIterations = 1000000;
        param.outputIntervalVTK = 10000;
        param.outputInterval = 10000;

        param.convergenceCriterium = 1.0e-10;

        param.L = 1.0;
        param.CFL = 0.1;

        param.verbose = false;
        param.fluxOutput = true;
        param.resOutput = false;

        // ========================================================================

        FluidParameter fluidParam;

        // ========== Parameters ==========
        int    ny = 3;//nyList[j];
        int    nx = 1;
        double Re = 40.0;//ReList[i];
        double u0 = 0.1;

        fluidParam.K = 1;
        fluidParam.nu = 0.1;
        fluidParam.R = 208.0;
        fluidParam.Force.x = 0.0;
        fluidParam.Force.y = 0.0;

        double TTop   = 273.15;
        double TBot   = 293.15;
        double lambda[] = { 1.0 / (2.0 * fluidParam.R * TTop), 1.0 / (2.0 * fluidParam.R * TBot) };
        double rho[]    = { 1.0, 1.0 * lambda[1] / lambda[0] };
        double U[] = { 0.0, 0.0 };
        double V[] = { 0.0, 0.0 };

        // ========================================================================

        GKSMesh* mesh = new GKSMesh(param, fluidParam);

        // Define Boundary Conditions
        //    -----------    -----------
        //    |    3    |    |    1    |
        //    | 0     2 |    |         |
        //    |    1    |    |    0    |
        //    -----------    -----------
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[0]);
        mesh->addBoundaryCondition(3, 0, 0, 0,  0.0, 0.0, 0.0, lambda[1]);

        Interface::setInterpolationOrder(1);

        // Generate Mesh
        mesh->generateRectMeshPeriodic(compressible, W, H, nx, ny);

        // Initialize Values
        mesh->initMeshLinear(rho, U, V, lambda);
        //mesh->initMeshConstant(1.0, 0.0, 0.0, 0.5*( lambda[0] + lambda[1] ) );

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
        //filename << "out/PoiseuilleFlowIncompressible/ConvergenceStudy1stOrder/" << "Re" << Re << "/" << ny;
        //mesh->writeVelocityProfile(            ( filename.str() + "/VelocityProfile.dat" )          , 1.0);
        //mesh->writePressureGradientProfile(    ( filename.str() + "/PressureGradientProfile.dat" )  , 1.0);
        //mesh->writeConvergenceHistory(         ( filename.str() + "/ConvergenceHistory.dat" )            );
        //mesh->writeOverviewFile(               ( filename.str() + "/OverviewFile.dat" )                  );
        //mesh->writeVTKFile(                    ( filename.str() + "/ResultFields.vtk" )                  );
    
        //ostringstream filename;
        //filename << "out/DrivenCavity/Re" << Re << "/" << nx;
        ////mesh->writeVelocityProfile(( filename.str() + "/VelocityProfile.dat" ), 0.5);
        //mesh->writeVelocityU(          (filename.str() + "/VelocityU.dat"          ));
        //mesh->writeVelocityV(          (filename.str() + "/VelocityV.dat"          ));
        //mesh->writeConvergenceHistory(( filename.str() + "/ConvergenceHistory.dat" ));
        //mesh->writeOverviewFile(      ( filename.str() + "/OverviewFile.dat" ));


        mesh->writeConvergenceHistory("out/ConvergenceHistory.dat");
        mesh->writeOverviewFile("out/OverviewFile.dat");
        mesh->writePressureGradientProfile("out/PressureGradientProfile.dat", 0.5);
        mesh->writeVelocityProfile("out/VelocityProfile.dat", 0.5);
        //mesh->writeTimeSteps("out/TimeSteps.dat");

        /*Cell* Cell1 = new Cell(compressible, -0.5, 0.5, 1.0, 1.0, NULL, fluidParam);
        Cell* Cell2 = new Cell(compressible,  0.5, 0.5, 1.0, 1.0, NULL, fluidParam);
        Cell* Cell3 = new Cell(compressible,  1.5, 0.5, 1.0, 1.0, NULL, fluidParam);
        Cell* Cell4 = new Cell(compressible,  2.5, 0.5, 1.0, 1.0, NULL, fluidParam);

        Interface* Interface1 = Interface::createInterface(compressible, Cell1, Cell2, float2( 0.0, 0.5 ), float2(1.0, 0.0), fluidParam, NULL);
        Interface* Interface2 = Interface::createInterface(compressible, Cell2, Cell3, float2( 0.0, 0.5 ), float2(1.0, 0.0), fluidParam, NULL);
        Interface* Interface3 = Interface::createInterface(compressible, Cell3, Cell4, float2( 0.0, 0.5 ), float2(1.0, 0.0), fluidParam, NULL);

        cout << Interface1->posCell->getOpposingCell(Interface1)->writeNodes() << endl;

        char a; cin >> a;*/
    }
}