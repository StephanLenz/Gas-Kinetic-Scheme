#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include "libgks/Reader/paramReader.h"

bool paramReader::read(string filename, Parameters & param, FluidParameter & fluidParam)
{
    paramReader::defaultValues(param, fluidParam);

    cout << "Start reading: " << filename << endl;
    ifstream file;
    file.open( filename.c_str() );

    if ( !file.is_open() ) {
        cout << " File cound not be opened.\n\nERROR!\n\n\n";
        return false;
    }

    string buffer;

    while( getline(file, buffer) )
    {
        stringstream bufferStream(buffer);

        string name, type;
        bufferStream >> name;

        if( name == "#" ) continue;

        if     ( name == "simulationName" )
        {
            bufferStream >> param.simulationName;
        }
        else if( name == "maxIter" )
        {
            bufferStream >> param.numberOfIterations;
        }
        else if( name == "maxTime" )
        {
            bufferStream >> param.maxTime;
        }
        else if( name == "output" )
        {
            bufferStream >> param.outputInterval;
        }
        else if( name == "outputVTK" )
        {
            bufferStream >> param.outputIntervalVTK;
        }
        else if( name == "convergenceCriterium" )
        {
            bufferStream >> param.convergenceCriterium.rho
                         >> param.convergenceCriterium.rhoU
                         >> param.convergenceCriterium.rhoV
                         >> param.convergenceCriterium.rhoE;
        }
        else if( name == "CFL" )
        {
            bufferStream >> param.CFL;
        }
        else if( name == "K" )
        {
            bufferStream >> fluidParam.K;
        }
        else if( name == "nu" )
        {
            bufferStream >> fluidParam.nu;
        }
        else if( name == "R" )
        {
            bufferStream >> fluidParam.R;
        }
        else if( name == "Pr" )
        {
            bufferStream >> fluidParam.Pr;
        }
        else if( name == "referencePrim" )
        {
            bufferStream >> fluidParam.referencePrim.rho
                         >> fluidParam.referencePrim.U
                         >> fluidParam.referencePrim.V
                         >> fluidParam.referencePrim.L;
        }
        else if( name == "Force" )
        {
            bufferStream >> fluidParam.Force.x
                         >> fluidParam.Force.y;
        }
        else if( name == "relForce" )
        {
            bufferStream >> fluidParam.relativeForce.x
                         >> fluidParam.relativeForce.y;
        }
        else
        {
            cout << "Error: Invalid Parameter " << type << endl;
            return false;
        }
    }

    file.close();

    cout << "Complete Parameters read!" << endl;

    return true;
}

void paramReader::defaultValues(Parameters& param, FluidParameter& fluidParam)
{
    param.numberOfIterations = 1000000;
    param.outputInterval     = 10000;
    param.outputIntervalVTK  = 10000;

    param.maxTime = 100.0;

    param.convergenceCriterium.rho  = 1.0e-10;
    param.convergenceCriterium.rhoU = 1.0e-10;
    param.convergenceCriterium.rhoV = 1.0e-10;
    param.convergenceCriterium.rhoE = 1.0e-10;

    param.CFL = 0.7;

    fluidParam.K  = 1.0;
    fluidParam.nu = 0.01;
    fluidParam.R  = 200;
    fluidParam.Pr = 1.0;

    fluidParam.referencePrim.rho = 1.0;
    fluidParam.referencePrim.U   = 0.0;
    fluidParam.referencePrim.V   = 0.0;
    fluidParam.referencePrim.L   = 0.0025;

    fluidParam.Force.x = 0.0;
    fluidParam.Force.y = 0.0;

    fluidParam.relativeForce.x = 0.0;
    fluidParam.relativeForce.y = 0.0;
}
