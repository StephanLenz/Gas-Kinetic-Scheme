
#include "libgks/Util/Types.h"
#include "libgks/Boundary/BoundaryCondition.h"
#include "libgks/Analyzer/FaceAnalyzer.h"
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

    static void defaultValues(Parameters& param, FluidParameter& fluidParam);
};

#endif
