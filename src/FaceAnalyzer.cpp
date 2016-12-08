#include "FaceAnalyzer.h"
#include "GKSSolver.h"
#include "mshReader.h"
#include <sstream>
#include <iostream>

using namespace std;

void FaceAnalyzer::addFace(idType id)
{
    this->Faces.push_back(id);
    this->Faces.shrink_to_fit();
}

DragAndLiftCalculator::DragAndLiftCalculator(double U, double D)
                      : drag(0.0), lift(0.0), U(U), D(D)
{
}

void DragAndLiftCalculator::analyze(GKSSolver& solver)
{
    double rhoSum = 0.0;
    double Fx = 0.0;
    double Fy = 0.0;

    for( int face = 0; face < this->Faces.size(); ++face )
    {
        ConservedVariable InterfaceFlux = this->computeTimeAveragedFlux( solver, this->Faces[face] );

        Fx -= InterfaceFlux.rhoU;
        Fy -= InterfaceFlux.rhoV;

        rhoSum += solver.getCellData( solver.getPosCell( face ) ).rho;
        rhoSum += solver.getCellData( solver.getNegCell( face ) ).rho;
    }

    double rho = rhoSum / ( 2.0 * this->Faces.size() );

    this->drag = ( 2.0 * Fx ) / ( rho * this->U * this->U * this->D );
    this->lift = ( 2.0 * Fy ) / ( rho * this->U * this->U * this->D );
}

void DragAndLiftCalculator::print()
{
    cout << "C_D = " << this->drag << "    C_L = " << this->lift << endl;
}

void DragAndLiftCalculator::write(string filename, double t)
{
    cout << "Wrinting file " << filename << " ... ";
	// open file stream
	ofstream file;
    file.precision(15);
	file.open(filename.c_str(), fstream::app);

	if (!file.is_open()) {
		cout << " File cound not be opened.\n\nERROR!\n\n\n";
		return;
	}

    file << t << ' ' << this->drag << ' ' << this->lift << endl;

    file.close();

    cout << "done!" << endl;
}

ConservedVariable DragAndLiftCalculator::computeTimeAveragedFlux(GKSSolver & solver, idType face)
{
    const int NUMBER_OF_MOMENTS = 7;

    PrimitiveVariable prim;
    ConservedVariable normalGradCons;
    ConservedVariable tangentialGradCons;
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
    prim = solver.reconstructPrimPiecewiseConstant(face);
    // ========================================================================

    // ========================================================================
    //          compute spacial gradients of the conservative varibles
    // ========================================================================
    //normalGradCons = differentiateConsNormal(id, prim.rho);
    solver.computeInterfaceGradient( face, prim.rho, normalGradCons, tangentialGradCons );
    // ========================================================================

    // ========================================================================
    //          Momentum Transformation in local coordinate system
    // ========================================================================
    solver.global2local(face, prim);
    solver.global2local(face, normalGradCons);
    solver.global2local(face, tangentialGradCons);
    // ========================================================================

    // ========================================================================
    //          compute spacial micro slopes
    //              a = a1 + a2 u + a3 v + 0.5 a4 (u^2 + v^2 + xi^2)
    // ========================================================================
    solver.computeMicroSlope(prim, normalGradCons, a);
    solver.computeMicroSlope(prim, tangentialGradCons, b);
    // ========================================================================

    // ========================================================================
    //          comoute moments of the equilibrium distribution
    // ========================================================================
    solver.computeMoments(prim, MomentU, MomentV, MomentXi, NUMBER_OF_MOMENTS);
    // ========================================================================

    // ========================================================================
    //          compute time derivative and temporal micro slopes
    //              A = A1 + A2 u + A3 v + 0.5 A4 (u^2 + v^2 + xi^2)
    // ========================================================================
    timeGrad = solver.computeTimeDerivative(MomentU, MomentV, MomentXi, a, b);

    solver.computeMicroSlope(prim, timeGrad, A);
    // ========================================================================

    // ========================================================================
    // Relaxation time as in the Rayleigh-Bernard-Paper (Xu, Lui, 1999)
    // ========================================================================
    double tau = 2.0*prim.L * solver.getFluidParam().nu;
    // ========================================================================

    // ========================================================================
    //          compute time integration Coefficients
    // ========================================================================
    double dt = solver.getDt();
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };
    // ========================================================================

    // ========================================================================
    //          compute mass and momentum fluxes
    // ========================================================================
    InterfaceFlux = solver.assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients, prim, solver.getInterfaceArea(face), tau);
    // ========================================================================



    // ========================================================================
    //          transform momentum Flux components back to global system
    // ========================================================================
    solver.local2global( face, InterfaceFlux );
    // ========================================================================
    InterfaceFlux.rho  /= dt;
    InterfaceFlux.rhoU /= dt;
    InterfaceFlux.rhoV /= dt;
    InterfaceFlux.rhoE /= dt;

    return InterfaceFlux;
}
