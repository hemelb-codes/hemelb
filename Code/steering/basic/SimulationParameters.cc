#include "io/XdrMemWriter.h"

#include "vis/Control.h"
#include "steering/basic/SimulationParameters.h"

namespace hemelb
{
  namespace steering
  {

    SimulationParameters::SimulationParameters() :
      paramWriter(params, paramsSizeB)
    {
      // C'tor initialises to the following defaults.

      pressureMin = 0.001;
      pressureMax = 1.0;
      velocityMin = 0.0;
      velocityMax = 1.0;
      stressMin = 0.0;
      stressMax = 1.0;
      timeStep = 0;
      time = 0.0;
      cycle = 0;
      nInlets = 3;
      mousePressure = -1.0;
      mouseStress = -1.0;
    }

    SimulationParameters::~SimulationParameters()
    {
    }

    char* SimulationParameters::pack()
    {
      paramWriter << pressureMin;
      paramWriter << pressureMax;

      paramWriter << velocityMin;
      paramWriter << velocityMax;

      paramWriter << stressMin;
      paramWriter << stressMax;

      paramWriter << timeStep;

      paramWriter << time;

      paramWriter << cycle;
      paramWriter << nInlets;

      paramWriter << mousePressure;
      paramWriter << mouseStress;

      return params;
    }

  }
}
