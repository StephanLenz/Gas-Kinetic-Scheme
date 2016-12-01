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

DragAndLiftCalculator::DragAndLiftCalculator()
{
}

void DragAndLiftCalculator::analyze(GKSSolver& solver)
{
}

void DragAndLiftCalculator::print()
{
    cout << "F_D = " << this->drag << "F_L = " << this->lift << endl;
}
