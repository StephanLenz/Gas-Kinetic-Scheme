// ============================================================================
//
//                      Compressible Thermal GKS
//
//      Developed by Stephan Lenz (stephan.lenz@tu-bs.de)
//
// ============================================================================
//
//      GKSSolver.h
//
//      Function:
//          Generation and Storage of mesh
//          Control of simulation
//          Data Analysis
//          File Output
//
// ============================================================================

#ifndef OUTPUTWRITER_H
#define OUTPUTWRITER_H

#include "GKSSolver.h"
#include "FaceAnalyzer.h"
#include "Types.h"
#include <fstream>

class outputWriter
{
private:
    outputWriter(){};
    outputWriter(outputWriter& orig){};

public:
    static void writeCellVTK(string filename, GKSSolver& solver);
    static void writeFaceVTK(string filename, GKSSolver& solver);

    static void initFile(string filename);
    static bool open(ofstream& file, string filename);

    static void writeOverview(string filename, GKSSolver& solver);

    static void writeConvergenceHistory(string filename, GKSSolver& solver);
};

#endif