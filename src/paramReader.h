
#include "Types.h"
#include "BoundaryCondition.h"
#include "FaceAnalyzer.h"
#include <vector>
#include <array>
#include <fstream>

#ifndef PARAMREADER_H
#define PARAMREADER_H

using namespace std;

class paramReader
{
public:
    static bool read(string filename, Parameters& param, FluidParameter& fluidParam);

    static void default(Parameters& param, FluidParameter& fluidParam);
};

#endif
