
#include "GKSSolver.h"
#include "Cell.h"
#include "Interface.h"
#include "BoundaryCondition.h"
#include "InterfaceBC.h"
#include "Types.h"
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

GKSSolver::GKSSolver() : dt(0.0), time(0.0), iter(0), computationTime(0)
{
}

GKSSolver::GKSSolver(Parameters param, FluidParameter fluidParam)
         : GKSSolver()
{
    this->fluidParam = fluidParam;
    this->param = param;
}

GKSSolver::~GKSSolver()
{
}

// ============================================================================
//      This method perform the iteration and controls the solution process. 
//      It also call the file output and checks the convergence.
// ============================================================================
void GKSSolver::iterate()
{
    this->iter = 0;
    this->time = 0.0;

    this->timeList.push_back(this->time);

    chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

    // ========================================================================
    // ========================================================================
    //              Time Loop
    // ========================================================================
    // ========================================================================
    while (this->iter < this->param.numberOfIterations)
    {
        // ====================================================================
        this->iter++;
        // ====================================================================

        // ====================================================================
        //          Perform one timestep
        // ====================================================================
        this->timeStep();
        // ====================================================================

        this->time += this->dt;
        
        // ====================================================================
        if ( this->iter % this->param.outputInterval == 0 )
        {
            this->dtList.push_back(this->dt);
            this->timeList.push_back(this->time);

            ConservedVariable residual = this->getL2GlobalResidual();

            cout << "t = " << this->time << "  \t|  timestep: " << this->iter << "  \t|  U_max: " << this->getMaxVelocity() << endl;
            cout << "r_rho = "  << residual.rho  << "\t ";
            cout << "r_rhoU = " << residual.rhoU << "\t ";
            cout << "r_rhoV = " << residual.rhoV << "\t ";
            cout << "r_rhoE = " << residual.rhoE << "\t ";
            cout << endl;

            this->convergenceHistory.push_back(residual);

            if ( this->isConverged(residual) )
            {
                cout << endl << " ========== Simulation converged! ==========" << endl;
                cout << "Remaining residual change less than " << this->param.convergenceCriterium << endl;
                cout << "Timesteps: " << this->iter << endl;
                cout << "Time: " << this->time << endl;

                //this->writeOutputFiles();

                this->writeVTK("out/SolverResult.vtk");

                break;
            }
        }
        // ====================================================================
        if ( this->iter % this->param.outputIntervalVTK == 0 )
        {
            this->writeVTK( string("out/Solver_") + to_string(this->iter) + string(".vtk") );
        }
    }
    // ========================================================================
    // ========================================================================
    //              End of Time Loop
    // ========================================================================
    // ========================================================================
    
    chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();
    this->computationTime = chrono::duration_cast<chrono::seconds>( endTime - startTime ).count();

    cout << "Time to Solution: " << this->computationTime << " s" << endl;
}

// ============================================================================
//      This method performes one time step of the GKS.
//      All computation is taking place with in this method.
// ============================================================================
void GKSSolver::timeStep()
{
    this->computeGlobalTimestep();
    
    #pragma omp parallel for
    for ( int id = 0; id < numberOfCells; ++id )
    {
        if ( !isGhostCell(id) ) applyForcing(id);
    }
    
    for( BoundaryCondition& BC : BoundaryConditionList )
    {
        BC.setGhostCells(*this);
    }

    //#pragma omp parallel for
    for ( int id = 0; id < numberOfInterfaces; ++id )
    {
        computeFlux(id);  
    }

    int breakPoint = 0;

    #pragma omp parallel for
    for ( int id = 0; id < numberOfCells; ++id )
    {
        if ( !isGhostCell(id) ) updateCell(id);   
    }

    breakPoint = 0;
}

bool GKSSolver::isConverged(ConservedVariable residual)
{
    bool flag = true;

    flag = flag && ( residual.rho  < this->param.convergenceCriterium[0] );
    flag = flag && ( residual.rhoU < this->param.convergenceCriterium[1] );
    flag = flag && ( residual.rhoV < this->param.convergenceCriterium[2] );
    flag = flag && ( residual.rhoE < this->param.convergenceCriterium[3] );

    return flag;
}

// ============================================================================
//      This method computes the global CFL time step by finding the minimum
//      of all local time steps and multiplying it with the CFL number.
// ============================================================================
void GKSSolver::computeGlobalTimestep()
{
    this->dt = 1.0e99;
    for (int id = 0; id < this->numberOfCells; ++id)
    {
        if ( ! isGhostCell(id) )
        {
            PrimitiveVariable prim = this->cons2prim( this->getCellData(id) );

            double U_max = sqrt(prim.U*prim.U + prim.V*prim.V);
            double c_s   = sqrt( 1.0 / ( 2.0*prim.L ) );           // c_s = sqrt(RT) = c_s = sqrt(1/2lambda)
            double minDx = this->getCellMinDx(id);

            double localTimestep = minDx / ( U_max + c_s + 2.0*this->fluidParam.nu / minDx );

            this->dt = min( localTimestep, this->dt );
        }
    }
    this->dt *= this->param.CFL;
}

// ============================================================================
//      This method applies the forcing to all cells
// ============================================================================
void GKSSolver::applyForcing(idType id)
{
    this->storeDataOld(id);
    
    // ========================================================================
    //                  Compute Primitive Variables (esp. Temp before forcing)
    // ========================================================================
    //PrimitiveVariable  prim = cons2prim( this->getCellData(id) );
    ConservedVariable cons = this->getCellData(id);
    PrimitiveVariable prim = cons2prim(cons);

    // ========================================================================
    //                  Apply Forcing to momentum components
    // ========================================================================
    cons.rhoU += dt * this->fluidParam.Force.x * cons.rho;
    cons.rhoV += dt * this->fluidParam.Force.y * cons.rho;
    // ========================================================================

    // ========================================================================
    //                  compute new total Energy with increased momentum
    // ========================================================================
    cons.rhoE = prim.rho * (this->fluidParam.K + 2.0) / (4.0*prim.L)
              + 0.5 * ( cons.rhoU*cons.rhoU + cons.rhoV*cons.rhoV ) / prim.rho;
    // ========================================================================

    this->setData(id, cons);
}

// ============================================================================
//      This method computes the fluxes over all Interfaces
// ============================================================================
void GKSSolver::computeFlux(const idType id)
{
    const int NUMBER_OF_MOMENTS = 7;

    PrimitiveVariable prim;
    ConservedVariable normalGradCons;
    ConservedVariable timeGrad;

    ConservedVariable InterfaceFlux;

    double a[4];
    double b[4] = {0.0, 0.0, 0.0, 0.0};
    double A[4];

    double MomentU[NUMBER_OF_MOMENTS];
    double MomentV[NUMBER_OF_MOMENTS];
    double MomentXi[NUMBER_OF_MOMENTS];

    // ========================================================================
    //          interpolated primitive variables at the interface
    // ========================================================================
    prim = reconstructPrimPiecewiseConstant(id);
    // ========================================================================
    
    // ========================================================================
    //          compute spacial gradients of the conservative varibles
    // ========================================================================
    normalGradCons = differentiateConsNormal(id, prim.rho);
    // ========================================================================
    
    // ========================================================================
    //          Momentum Transformation in local coordinate system
    // ========================================================================
    global2local(id, prim);
    global2local(id, normalGradCons);
    // ========================================================================
    
    // ========================================================================
    //          compute spacial micro slopes
    //              a = a1 + a2 u + a3 v + 0.5 a4 (u^2 + v^2 + xi^2)
    // ========================================================================
    this->computeMicroSlope(prim, normalGradCons, a);
    // ========================================================================
    
    // ========================================================================
    //          comoute moments of the equilibrium distribution
    // ========================================================================
    this->computeMoments(prim, MomentU, MomentV, MomentXi, NUMBER_OF_MOMENTS);
    // ========================================================================

    // ========================================================================
    //          compute time derivative and temporal micro slopes
    //              A = A1 + A2 u + A3 v + 0.5 A4 (u^2 + v^2 + xi^2)
    // ========================================================================
    timeGrad = computeTimeDerivative(MomentU, MomentV, MomentXi, a, b);

    this->computeMicroSlope(prim, timeGrad, A);
    // ========================================================================

    // ========================================================================
    // Relaxation time as in the Rayleigh-Bernard-Paper (Xu, Lui, 1999)
    // ========================================================================
    double tau = 2.0*prim.L * this->fluidParam.nu;
    // ========================================================================
    
    // ========================================================================
    //          compute time integration Coefficients
    // ========================================================================
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };
    // ========================================================================

    // ========================================================================
    //          compute mass and momentum fluxes
    // ========================================================================
    InterfaceFlux = assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients, prim, getInterfaceArea(id), tau);
    // ========================================================================

    // ========================================================================
    //          transform momentum Flux components back to global system
    // ========================================================================
    local2global( id, InterfaceFlux );
    // ========================================================================
    //#pragma omp ordered
    //#pragma omp critical
    this->applyFlux(id, InterfaceFlux);
}

PrimitiveVariable GKSSolver::reconstructPrimPiecewiseConstant(const idType id)
{
    PrimitiveVariable posPrim = cons2prim( this->getCellData( this->getPosCell(id) ) );
    PrimitiveVariable negPrim = cons2prim( this->getCellData( this->getNegCell(id) ) );

    PrimitiveVariable midPrim;

    midPrim.rho  = 0.5 * ( posPrim.rho + negPrim.rho );
    midPrim.U    = 0.5 * ( posPrim.U   + negPrim.U   );
    midPrim.V    = 0.5 * ( posPrim.V   + negPrim.V   );
    midPrim.L    = 0.5 * ( posPrim.L   + negPrim.L   );

    return midPrim;
}

ConservedVariable GKSSolver::differentiateConsNormal(const idType id, double rho)
{
    ConservedVariable posCons = this->getCellData( this->getPosCell(id) );
    ConservedVariable negCons = this->getCellData( this->getNegCell(id) );

    ConservedVariable gradConsNormal;
    gradConsNormal.rho  = ( posCons.rho  - negCons.rho  ) / ( this->getInterfaceDistance(id) * rho );
    gradConsNormal.rhoU = ( posCons.rhoU - negCons.rhoU ) / ( this->getInterfaceDistance(id) * rho );
    gradConsNormal.rhoV = ( posCons.rhoV - negCons.rhoV ) / ( this->getInterfaceDistance(id) * rho );
    gradConsNormal.rhoE = ( posCons.rhoE - negCons.rhoE ) / ( this->getInterfaceDistance(id) * rho );

    return gradConsNormal;
}

void GKSSolver::computeMicroSlope(const PrimitiveVariable& prim, const ConservedVariable& macroSlope, double * microSlope)
{
    // this method computes the micro slopes from the slopes of the conservative variables
    // the resulting microslopes contain the density, since they are computed from the slopes
    // of the conservative variables, which are rho, rhoU, rhoV and rhoE

    double A, B, C, E;

    // ========================================================================
    // this is 2 times the total energy density 2 E = 2 rhoE / rho
    E = prim.U * prim.U + prim.V * prim.V + ( this->fluidParam.K + 2.0 ) / ( 2.0*prim.L );
    // ========================================================================

    // ========================================================================
    // the product rule of derivations is used here!
    A = 2.0*macroSlope.rhoE -      E  * macroSlope.rho;    // = 2 rho dE/dx
    B =     macroSlope.rhoU - prim.U  * macroSlope.rho;    // =   rho dU/dx
    C =     macroSlope.rhoV - prim.V  * macroSlope.rho;    // =   rho dV/dx
    // ========================================================================

    // compute micro slopes of primitive variables from macro slopes of conservative variables
    microSlope[3] = ( 4.0 * prim.L * prim.L ) / ( this->fluidParam.K + 2.0 )
                  * ( A - 2.0*prim.U * B - 2.0*prim.V * C );

    microSlope[2] = 2.0 * prim.L * C - prim.V * microSlope[3];

    microSlope[1] = 2.0 * prim.L * B - prim.U * microSlope[3];

    microSlope[0] = macroSlope.rho - prim.U * microSlope[1] - prim.V * microSlope[2] - 0.5 * E * microSlope[3];

}

void GKSSolver::computeMoments(const PrimitiveVariable & prim, double * MomentU, double * MomentV, double * MomentXi, const int numberMoments)
{
    //==================== U Moments ==========================================
    MomentU[0] = 1.0;
    MomentU[1] = prim.U;
    for (int i = 2; i < numberMoments; i++)
        MomentU[i] = prim.U * MomentU[i - 1] + ( double(i - 1) * MomentU[i - 2] )/( 2.0 * prim.L );

    //==================== V Moments ==========================================
    MomentV[0] = 1.0;
    MomentV[1] = prim.V;
    for (int i = 2; i < numberMoments; i++)
        MomentV[i] = prim.V * MomentV[i - 1] + ( double(i - 1) * MomentV[i - 2] )/( 2.0 * prim.L );

    //==================== Xi Moments =========================================
    MomentXi[0] = 1.0;
    MomentXi[1] = 0.0;
    MomentXi[2] = this->fluidParam.K / (2.0 * prim.L);
    MomentXi[3] = 0.0;
    MomentXi[4] = ( 2.0*this->fluidParam.K + this->fluidParam.K*this->fluidParam.K ) / (4.0 * prim.L * prim.L);
    MomentXi[5] = 0.0;
    MomentXi[6] = ( this->fluidParam.K + 4.0 ) / ( 2.0 * prim.L ) * MomentXi[4];
}

ConservedVariable GKSSolver::computeTimeDerivative(double * MomentU, double * MomentV, double * MomentXi, double * a, double * b)
{
    ConservedVariable timeGrad;

    // ========================================================================
    timeGrad.rho = a[0] * MomentU[1] * MomentV[0]
                 + a[1] * MomentU[2] * MomentV[0]
                 + a[2] * MomentU[1] * MomentV[1]
                 + a[3] * 0.5 * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                 + b[0] * MomentU[0] * MomentV[1]
                 + b[1] * MomentU[1] * MomentV[1]
                 + b[2] * MomentU[0] * MomentV[2]
                 + b[3] * 0.5 * ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] )
    // this part comes from the inclusion of the forcing into the flux computation
                 //+ 2.0 * prim[3] * ( MomentU[0]*MomentV[0] * prim[1] - MomentU[1]*MomentV[0] ) * this->fluidParam.Force.x
                 //+ 2.0 * prim[3] * ( MomentU[0]*MomentV[0] * prim[2] - MomentU[0]*MomentV[1] ) * this->fluidParam.Force.y
                 ;
    // ========================================================================
    
    // ========================================================================
    timeGrad.rhoU = a[0] * MomentU[2] * MomentV[0]
                  + a[1] * MomentU[3] * MomentV[0]
                  + a[2] * MomentU[2] * MomentV[1]
                  + a[3] * 0.5 * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                  + b[0] * MomentU[1] * MomentV[1]
                  + b[1] * MomentU[2] * MomentV[1]
                  + b[2] * MomentU[1] * MomentV[2]
                  + b[3] * 0.5 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
    // this part comes from the inclusion of the forcing into the flux computation
                  //+ 2.0 * prim[3] * ( MomentU[1]*MomentV[0] * prim[1] - MomentU[2]*MomentV[0] ) * this->fluidParam.Force.x
                  //+ 2.0 * prim[3] * ( MomentU[1]*MomentV[0] * prim[2] - MomentU[1]*MomentV[1] ) * this->fluidParam.Force.y
                  ;
    // ========================================================================
    
    // ========================================================================
    timeGrad.rhoV = a[0] * MomentU[1] * MomentV[1]
                  + a[1] * MomentU[2] * MomentV[1]
                  + a[2] * MomentU[1] * MomentV[2]
                  + a[3] * 0.5 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                  + b[0] * MomentU[0] * MomentV[2]
                  + b[1] * MomentU[1] * MomentV[2]
                  + b[2] * MomentU[0] * MomentV[3]
                  + b[3] * 0.5 * ( MomentU[2] * MomentV[2] + MomentU[0] * MomentV[4] + MomentU[0] * MomentV[2] * MomentXi[2] )
    // this part comes from the inclusion of the forcing into the flux computation
                  //+ 2.0 * prim[3] * ( MomentU[0]*MomentV[1] * prim[1] - MomentU[1]*MomentV[1] ) * this->fluidParam.Force.x
                  //+ 2.0 * prim[3] * ( MomentU[0]*MomentV[1] * prim[2] - MomentU[0]*MomentV[2] ) * this->fluidParam.Force.y
                  ;
    // ========================================================================
    
    // ========================================================================
    timeGrad.rhoE = a[0] * 0.50 * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                  + a[1] * 0.50 * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                  + a[2] * 0.50 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                  + a[3] * 0.25 * ( MomentU[5] + MomentU[1]* ( MomentV[4] + MomentXi[4] )
                                  + 2.0 * MomentU[3] * MomentV[2]
                                  + 2.0 * MomentU[3] * MomentXi[2]
                                  + 2.0 * MomentU[1] * MomentV[2] * MomentXi[2]
                                  )
                  + b[0] * 0.50 * ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] )
                  + b[1] * 0.50 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                  + b[2] * 0.50 * ( MomentU[2] * MomentV[2] + MomentU[0] * MomentV[4] + MomentU[0] * MomentV[2] * MomentXi[2] )
                  + b[3] * 0.25 * ( MomentV[5] + MomentV[1] * ( MomentU[4] + MomentXi[4] )
                                  + 2.0 * MomentU[2] * MomentV[3]
                                  + 2.0 * MomentU[2] * MomentV[1] * MomentXi[2]
                                  + 2.0 * MomentV[3] * MomentXi[2]
                                  )
    // this part comes from the inclusion of the forcing into the flux computation
                  //+ prim[3] * ( ( MomentU[2] + MomentV[2] + MomentXi[2] ) * prim[1] 
                  //            - ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                  //            ) * this->fluidParam.Force.x
                  //+ prim[3] * ( ( MomentU[2] + MomentV[2] + MomentXi[2] ) * prim[2] 
                  //            - ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] ) 
                  //            ) * this->fluidParam.Force.y
                  ;
    // ========================================================================

    timeGrad.rho  *= -1.0;
    timeGrad.rhoU *= -1.0;
    timeGrad.rhoV *= -1.0;
    timeGrad.rhoE *= -1.0;

    return timeGrad;
}

ConservedVariable GKSSolver::assembleFlux(double * MomentU, double * MomentV, double * MomentXi, double * a, double * b, double * A, double * timeCoefficients, PrimitiveVariable prim, double area, double tau)
{
    ConservedVariable Flux_1;
    ConservedVariable Flux_2;
    ConservedVariable Flux_3;

    ConservedVariable Flux;
    
    // ================================================================================================================================================
    Flux_1.rho  = MomentU[0+1] * MomentV[0];
    Flux_1.rhoU = MomentU[1+1] * MomentV[0];
    Flux_1.rhoV = MomentU[0+1] * MomentV[1];
    Flux_1.rhoE = 0.5 * ( MomentU[2+1] * MomentV[0]
                        + MomentU[0+1] * MomentV[2] 
                        + MomentU[0+1] * MomentV[0] * MomentXi[2] );
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================

    // ================================================================================================================================================
    Flux_2.rho = ( a[0] * MomentU[1+1] * MomentV[0]
                 + a[1] * MomentU[2+1] * MomentV[0]
                 + a[2] * MomentU[1+1] * MomentV[1]
                 + a[3] * 0.5 * ( MomentU[3+1] * MomentV[0]
                                + MomentU[1+1] * MomentV[2] 
                                + MomentU[1+1] * MomentV[0] * MomentXi[2] )
                 + b[0] * MomentU[0+1] * MomentV[1]
                 + b[1] * MomentU[1+1] * MomentV[1]
                 + b[2] * MomentU[0+1] * MomentV[2]
                 + b[3] * 0.5 * ( MomentU[2+1] * MomentV[1] 
                                + MomentU[0+1] * MomentV[3] 
                                + MomentU[0+1] * MomentV[1] * MomentXi[2] )
                 )
    // this part comes from the inclusion of the forcing into the flux computation
               //+ 2.0 * prim[3] * ( MomentU[0+1]*MomentV[0] * prim[1] - MomentU[1+1]*MomentV[0] ) * this->fluidParam.Force.x
               //+ 2.0 * prim[3] * ( MomentU[0+1]*MomentV[0] * prim[2] - MomentU[0+1]*MomentV[1] ) * this->fluidParam.Force.y
               ;
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    Flux_2.rhoU = ( a[0] * MomentU[2+1] * MomentV[0]
                  + a[1] * MomentU[3+1] * MomentV[0]
                  + a[2] * MomentU[2+1] * MomentV[1]
                  + a[3] * 0.5 * ( MomentU[4+1] * MomentV[0] 
                                 + MomentU[2+1] * MomentV[2] 
                                 + MomentU[2+1] * MomentV[0] * MomentXi[2] )
                  + b[0] * MomentU[1+1] * MomentV[1]
                  + b[1] * MomentU[2+1] * MomentV[1]
                  + b[2] * MomentU[1+1] * MomentV[2]
                  + b[3] * 0.5 * ( MomentU[3+1] * MomentV[1] 
                                 + MomentU[1+1] * MomentV[3] 
                                 + MomentU[1+1] * MomentV[1] * MomentXi[2] )
                  )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ 2.0 * prim[3] * ( MomentU[1+1]*MomentV[0] * prim[1] - MomentU[2+1]*MomentV[0] ) * this->fluidParam.Force.x
                //+ 2.0 * prim[3] * ( MomentU[1+1]*MomentV[0] * prim[2] - MomentU[1+1]*MomentV[1] ) * this->fluidParam.Force.y
                ;
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    Flux_2.rhoV = ( a[0] * MomentU[1+1] * MomentV[1]
                  + a[1] * MomentU[2+1] * MomentV[1]
                  + a[2] * MomentU[1+1] * MomentV[2]
                  + a[3] * 0.5 * ( MomentU[3+1] * MomentV[1]
                                 + MomentU[1+1] * MomentV[3]
                                 + MomentU[1+1] * MomentV[1] * MomentXi[2] )
                  + b[0] * MomentU[0+1] * MomentV[2]
                  + b[1] * MomentU[1+1] * MomentV[2]
                  + b[2] * MomentU[0+1] * MomentV[3]
                  + b[3] * 0.5 * ( MomentU[2+1] * MomentV[2]
                                 + MomentU[0+1] * MomentV[4]
                                 + MomentU[0+1] * MomentV[2] * MomentXi[2] )
                  )
    // this part comes from the inclusion of the forcing into the flux computation
                //+ 2.0 * prim[3] * ( MomentU[0+1]*MomentV[1] * prim[1] - MomentU[1+1]*MomentV[1] ) * this->fluidParam.Force.x
                //+ 2.0 * prim[3] * ( MomentU[0+1]*MomentV[1] * prim[2] - MomentU[0+1]*MomentV[2] ) * this->fluidParam.Force.y
                ;
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    Flux_2.rhoE = 0.5 * ( a[0] * ( MomentU[3+1] * MomentV[0] 
                                 + MomentU[1+1] * MomentV[2]
                                 + MomentU[1+1] * MomentV[0] * MomentXi[2] )
                        + a[1] * ( MomentU[4+1] * MomentV[0]
                                 + MomentU[2+1] * MomentV[2]
                                 + MomentU[2+1] * MomentV[0] * MomentXi[2] )
                        + a[2] * ( MomentU[3+1] * MomentV[1]
                                 + MomentU[1+1] * MomentV[3]
                                 + MomentU[1+1] * MomentV[1] * MomentXi[2] )
                        + a[3] * ( 0.5 * ( MomentU[5+1] * MomentV[0]
                                         + MomentU[1+1] * MomentV[4]
                                         + MomentU[1+1] * MomentV[0] * MomentXi[4] )
                                 +       ( MomentU[3+1] * MomentV[2]
                                         + MomentU[3+1] * MomentV[0] * MomentXi[2]
                                         + MomentU[1+1] * MomentV[2] * MomentXi[2] ) 
                                 )
                        + b[0] * ( MomentU[2+1] * MomentV[1] 
                                 + MomentU[0+1] * MomentV[3]
                                 + MomentU[0+1] * MomentV[1] * MomentXi[2] )
                        + b[1] * ( MomentU[3+1] * MomentV[1]
                                 + MomentU[1+1] * MomentV[3]
                                 + MomentU[1+1] * MomentV[1] * MomentXi[2] )
                        + b[2] * ( MomentU[2+1] * MomentV[2]
                                 + MomentU[0+1] * MomentV[4]
                                 + MomentU[0+1] * MomentV[2] * MomentXi[2] )
                        + b[3] * ( 0.5 * ( MomentU[4+1] * MomentV[1] 
                                         + MomentU[0+1] * MomentV[5]
                                         + MomentU[0+1] * MomentV[1] * MomentXi[4] )
                                 +       ( MomentU[2+1] * MomentV[3]
                                         + MomentU[2+1] * MomentV[1] * MomentXi[2]
                                         + MomentU[0+1] * MomentV[3] * MomentXi[2] )
                                 )
                        )
    // this part comes from the inclusion of the forcing into the flux computation
                      //+ prim[3] * ( ( MomentU[2+1] * MomentV[0] + MomentU[0+1] * MomentV[2] + MomentU[0+1] * MomentV[0] * MomentXi[2] ) * prim[1] 
                      //            - ( MomentU[3+1] * MomentV[0] + MomentU[1+1] * MomentV[2] + MomentU[1+1] * MomentV[0] * MomentXi[2] )
                      //            ) * this->fluidParam.Force.x
                      //+ prim[3] * ( ( MomentU[2+1] * MomentV[0] + MomentU[0+1] * MomentV[2] + MomentU[0+1] * MomentV[0] * MomentXi[2] ) * prim[2] 
                      //            - ( MomentU[2+1] * MomentV[1] + MomentU[0+1] * MomentV[3] + MomentU[0+1] * MomentV[1] * MomentXi[2] ) 
                      //            ) * this->fluidParam.Force.y
                      ;
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================

    // ================================================================================================================================================
    Flux_3.rho = ( A[0] * MomentU[0+1] * MomentV[0]
                 + A[1] * MomentU[1+1] * MomentV[0]
                 + A[2] * MomentU[0+1] * MomentV[1]
                 + A[3] * 0.5 * ( MomentU[2+1]*MomentV[0]
                                + MomentU[0+1]*MomentV[2]
                                + MomentU[0+1]*MomentV[0]*MomentXi[2] )
                 );
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    Flux_3.rhoU = ( A[0] * MomentU[1+1] * MomentV[0]
                  + A[1] * MomentU[2+1] * MomentV[0]
                  + A[2] * MomentU[1+1] * MomentV[1]
                  + A[3] * 0.5 * ( MomentU[3+1]*MomentV[0]
                                 + MomentU[1+1]*MomentV[2]
                                 + MomentU[1+1]*MomentV[0]*MomentXi[2] )
                  );
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    Flux_3.rhoV = ( A[0] * MomentU[0+1] * MomentV[1]
                  + A[1] * MomentU[1+1] * MomentV[1]
                  + A[2] * MomentU[0+1] * MomentV[2]
                  + A[3] * 0.5 * ( MomentU[2+1]*MomentV[1]
                                 + MomentU[0+1]*MomentV[3]
                                 + MomentU[0+1]*MomentV[1]*MomentXi[2] )
                  );
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    Flux_3.rhoE = 0.5 * ( A[0] * ( MomentU[2+1] * MomentV[0]
                                 + MomentU[0+1] * MomentV[2]
                                 + MomentU[0+1] * MomentV[0] * MomentXi[2] )
                        + A[1] * ( MomentU[3+1] * MomentV[0]
                                 + MomentU[1+1] * MomentV[2]
                                 + MomentU[1+1] * MomentV[0] * MomentXi[2] )
                        + A[2] * ( MomentU[2+1] * MomentV[1]
                                 + MomentU[0+1] * MomentV[3]
                                 + MomentU[0+1] * MomentV[1] * MomentXi[2] )
                        + A[3] * ( 0.5 * ( MomentU[4+1] * MomentV[0]
                                         + MomentU[0+1] * MomentV[4]
                                         + MomentU[0+1] * MomentV[0] * MomentXi[4] )
                                 +       ( MomentU[2+1] * MomentV[2]
                                         + MomentU[2+1] * MomentV[0] * MomentXi[2]
                                         + MomentU[0+1] * MomentV[2] * MomentXi[2] ) 
                                 )
                        );
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    // ================================================================================================================================================
    // ================================================================================================================================================

    // ================================================================================================================================================
    Flux.rho  = ( timeCoefficients[0] * Flux_1.rho  + timeCoefficients[1] * Flux_2.rho  + timeCoefficients[2] * Flux_3.rho  ) * area * prim.rho;
    Flux.rhoU = ( timeCoefficients[0] * Flux_1.rhoU + timeCoefficients[1] * Flux_2.rhoU + timeCoefficients[2] * Flux_3.rhoU ) * area * prim.rho;
    Flux.rhoV = ( timeCoefficients[0] * Flux_1.rhoV + timeCoefficients[1] * Flux_2.rhoV + timeCoefficients[2] * Flux_3.rhoV ) * area * prim.rho;
    Flux.rhoE = ( timeCoefficients[0] * Flux_1.rhoE + timeCoefficients[1] * Flux_2.rhoE + timeCoefficients[2] * Flux_3.rhoE ) * area * prim.rho;
    // ================================================================================================================================================

    // ================================================================================================================================================
    //                              Prandl number correction
    // ================================================================================================================================================
    double q = timeCoefficients[1] * ( ( Flux_2.rhoE + Flux_3.rhoE ) - prim.U * ( Flux_2.rhoU + Flux_3.rhoU ) - prim.V * ( Flux_2.rhoV + Flux_3.rhoV ) ) * area * prim.rho;

    Flux.rhoE += ( 1.0/this->fluidParam.Pr - 1.0 ) * q;
    // ================================================================================================================================================
    
    // ================================================================================================================================================
    //                              Pressure Conditioning
    // ================================================================================================================================================
    //double pReference = this->fluidParam.rhoReference / ( 2.0 * this->fluidParam.lambdaReference );

    //Flux.rhoU -= pReference * this->dt * area;

    return Flux;
}

// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================
//
//                          Data Analysis
//
// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

// ============================================================================
//      This method returns the relative L2 norm of the residual change per time step.
// ============================================================================
ConservedVariable GKSSolver::getL2GlobalResidual()
{
    ConservedVariable residual;
    ConservedVariable residualSquare;
    residualSquare.rho  = 0.0;
    residualSquare.rhoU = 0.0;
    residualSquare.rhoV = 0.0;
    residualSquare.rhoE = 0.0;

    ConservedVariable consSquare;
    consSquare.rho  = 0.0;
    consSquare.rhoU = 0.0;
    consSquare.rhoV = 0.0;
    consSquare.rhoE = 0.0;

    for ( int id = 0; id < this->numberOfCells; ++id )
    {
        if ( ! isGhostCell(id) )
        {
            residualSquare.rho  +=  ( this->getCellDataOld(id).rho  - this->getCellData(id).rho  ) * ( this->getCellDataOld(id).rho  - this->getCellData(id).rho  );
            residualSquare.rhoU +=  ( this->getCellDataOld(id).rhoU - this->getCellData(id).rhoU ) * ( this->getCellDataOld(id).rhoU - this->getCellData(id).rhoU );
            residualSquare.rhoV +=  ( this->getCellDataOld(id).rhoV - this->getCellData(id).rhoV ) * ( this->getCellDataOld(id).rhoV - this->getCellData(id).rhoV );
            residualSquare.rhoE +=  ( this->getCellDataOld(id).rhoE - this->getCellData(id).rhoE ) * ( this->getCellDataOld(id).rhoE - this->getCellData(id).rhoE );

            consSquare.rho  +=  this->getCellData(id).rho  * this->getCellData(id).rho ;
            consSquare.rhoU +=  this->getCellData(id).rhoU * this->getCellData(id).rhoU;
            consSquare.rhoV +=  this->getCellData(id).rhoV * this->getCellData(id).rhoV;
            consSquare.rhoE +=  this->getCellData(id).rhoE * this->getCellData(id).rhoE;
        }
    }

    residual.rho  = sqrt( residualSquare.rho  ) / sqrt( consSquare.rho  );
    residual.rhoU = sqrt( residualSquare.rhoU ) / sqrt( consSquare.rhoU );
    residual.rhoV = sqrt( residualSquare.rhoV ) / sqrt( consSquare.rhoV );
    residual.rhoE = sqrt( residualSquare.rhoE ) / sqrt( consSquare.rhoE );

    return residual;
}

double GKSSolver::getMaxVelocity()
{
    double maxVelocity = 0.0;
    double localVelocity;
    for ( int id = 0; id < this->numberOfCells; ++id )
    {
        if ( ! isGhostCell(id) )
        {
            PrimitiveVariable prim = cons2prim( this->getCellData(id) );
            localVelocity = sqrt( prim.U * prim.U 
                                + prim.V * prim.V );
            maxVelocity = max(maxVelocity, localVelocity);
        }
    }

    return maxVelocity;
}

PrimitiveVariable GKSSolver::cons2prim(ConservedVariable& cons)
{
    PrimitiveVariable prim;
    prim.rho = cons.rho;
    prim.U   = cons.rhoU / cons.rho;
    prim.V   = cons.rhoV / cons.rho;
    prim.L = (this->fluidParam.K + 2.0)*cons.rho
           / ( 4.0 * ( cons.rhoE - 0.5*( cons.rhoU*cons.rhoU + cons.rhoV*cons.rhoV) / cons.rho ) );

    return prim;
}

ConservedVariable GKSSolver::prim2cons(PrimitiveVariable & prim)
{
    ConservedVariable cons;

    cons.rho  = prim.rho;
    cons.rhoU = prim.rho * prim.U; 
    cons.rhoV = prim.rho * prim.V;
    cons.rhoE = prim.rho * (this->fluidParam.K + 2.0) / (4.0*prim.L) + 0.5 * prim.rho * ( prim.U*prim.U + prim.V*prim.V );

    return cons;
}

void GKSSolver::global2local(const idType id, PrimitiveVariable & prim)
{
    // euclidian components in global coordinatesystem
    double u0 = prim.U;
    double v0 = prim.V;

    // transformation in local coordinatesystem
    // n = (n1,n2)
    // t = (-n2,n1)
    // vL = [n t]^T * v0
    prim.U =   getInterfaceNormal(id).x * u0 + getInterfaceNormal(id).y * v0;
    prim.V = - getInterfaceNormal(id).y * u0 + getInterfaceNormal(id).x * v0;
}

void GKSSolver::global2local(const idType id, ConservedVariable & cons)
{
    // euclidian components in global coordinatesystem
    double u0 = cons.rhoU;
    double v0 = cons.rhoV;

    // transformation in local coordinatesystem
    // n = (n1,n2)
    // t = (-n2,n1)
    // vL = [n t]^T * v0
    cons.rhoU =   getInterfaceNormal(id).x * u0 + getInterfaceNormal(id).y * v0;
    cons.rhoV = - getInterfaceNormal(id).y * u0 + getInterfaceNormal(id).x * v0;
}

void GKSSolver::local2global(const idType id, PrimitiveVariable & prim)
{
    // euclidian components in local coordinatesystem
    double un = prim.U;
    double vt = prim.V;

    // transformation in global coordinatesystem
    // n = (n1,n2)
    // t = (-n2,n1)
    // vL = [n t] * v0
    prim.U = getInterfaceNormal(id).x * un - getInterfaceNormal(id).y * vt;
    prim.V = getInterfaceNormal(id).y * un + getInterfaceNormal(id).x * vt;
}

void GKSSolver::local2global(const idType id, ConservedVariable & cons)
{
    // euclidian components in local coordinatesystem
    double un = cons.rhoU;
    double vt = cons.rhoV;

    // transformation in global coordinatesystem
    // n = (n1,n2)
    // t = (-n2,n1)
    // vL = [n t] * v0
    cons.rhoU = getInterfaceNormal(id).x * un - getInterfaceNormal(id).y * vt;
    cons.rhoV = getInterfaceNormal(id).y * un + getInterfaceNormal(id).x * vt;
}

void GKSSolver::setData(idType id, PrimitiveVariable prim)
{
    this->setData( id, prim2cons(prim) );
}

PrimitiveVariable GKSSolver::getPrim(idType id)
{
    return cons2prim( this->getCellData(id) );
}

void GKSSolver::writeVTK(string filename)
{
    cout << "Wrinting file " << filename << " ... ";
	// open file stream
	ofstream file;
    file.precision(15);
	file.open(filename.c_str());

	if (!file.is_open()) {
		cout << " File cound not be opened.\n\nERROR!\n\n\n";
		return;
	}
    
    // write VTK Header
    file << "# vtk DataFile Version 1.0\n";
    file << "by Stephan Lenz\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    
    file << "POINTS " << this->numberOfNodes << " double\n";
    for(int node = 0; node < this->numberOfNodes; ++node)
    {
        file << this->getNode(node).x << " " << this->getNode(node).y << " 0.0" << endl;
    }

    file << "CELLS " << numberOfCells << " " << 5 * numberOfCells << endl;
    for (int cell = 0; cell < this->numberOfCells; ++cell)
    {
        file << "4 " << this->getCell2Node(cell)[0] << " "
                     << this->getCell2Node(cell)[1] << " "
                     << this->getCell2Node(cell)[2] << " ";

        if( -1 == this->getCell2Node(cell)[3] )
            file << this->getCell2Node(cell)[2] << endl;
        else
            file << this->getCell2Node(cell)[3] << endl;
    }

    file << "CELL_TYPES " << numberOfCells << endl;
    for (int cell = 0; cell < this->numberOfCells; ++cell)
    {
        file << "9" << endl;
    }
    
    file << "CELL_DATA " << this->numberOfCells << endl;
    file << "FIELD Lable " << 5 << "\n";
    file << "rho 1 " << numberOfCells << " double\n";
    for (int cell = 0; cell < this->numberOfCells; ++cell)
    {
        file << this->getPrim(cell).rho << endl;
    }
    file << "U 1 " << numberOfCells << " double\n";
    for (int cell = 0; cell < this->numberOfCells; ++cell)
    {
        file << this->getPrim(cell).U << endl;
    }
    file << "V 1 " << numberOfCells << " double\n";
    for (int cell = 0; cell < this->numberOfCells; ++cell)
    {
        file << this->getPrim(cell).V << endl;
    }
    file << "T 1 " << numberOfCells << " double\n";
    for (int cell = 0; cell < this->numberOfCells; ++cell)
    {
        file << 1.0 / ( 2.0 * this->fluidParam.R * this->getPrim(cell).L ) << endl;
    }
    file << "Lambda 1 " << numberOfCells << " double\n";
    for (int cell = 0; cell < this->numberOfCells; ++cell)
    {
        file << this->getPrim(cell).L << endl;
    }

    file.close();

    cout << "done!" << endl;
}
