
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/stents/StentFileFlux.h"
#include <algorithm>
#include <fstream>
#include "log/Logger.h"
#include "util/fileutils.h"
#include "util/utilityFunctions.h"
#include "util/utilityStructs.h"
#include "configuration/SimConfig.h"
#include <cmath>
#include <algorithm>

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {
      StentFileFlux::StentFileFlux() :
          units(NULL)
      {
      }

      Stent* StentFileFlux::Clone() const
      {
        Stent* copy = new StentFileFlux(*this);
        return copy;
      }

      void StentFileFlux::CalculateTable(LatticeTimeStep totalTimeSteps, PhysicalTime timeStepLength)
      {
        // First read in values from file
        // Used to be complex code here to keep a vector unique, but this is just achieved by using a map.
        std::map<PhysicalTime, PhysicalSpeed> timeValuePairs;

        double timeTemp, valueTemp;

        util::check_file(fluxFilePath.c_str());
        std::ifstream datafile(fluxFilePath.c_str());
        log::Logger::Log<log::Debug, log::OnePerCore>("Reading stent values from file:");
        while (datafile.good())
        {
          datafile >> timeTemp >> valueTemp;
          log::Logger::Log<log::Trace, log::OnePerCore>("Time: %f Value: %f", timeTemp, valueTemp);
          timeValuePairs[timeTemp] = valueTemp;
        }

        datafile.close();
        // the default iterator for maps traverses in key order, so no sort is needed.

        std::vector<PhysicalTime> times(0);
        std::vector<PhysicalSpeed> values(0);

        // Must convert into vectors since LinearInterpolate works on a pair of vectors
        // Determine min and max pressure on the way
//        PhysicalPressure pMin = timeValuePairs.begin()->second;
//        PhysicalPressure pMax = timeValuePairs.begin()->second;
        for (std::map<PhysicalTime, PhysicalSpeed>::iterator entry = timeValuePairs.begin();
            entry != timeValuePairs.end(); entry++)
        {

          /* If the time value in the input file stretches BEYOND the end of the simulation, then insert an interpolated end value and exit the loop. */
          if(entry->first > totalTimeSteps*timeStepLength) {
  
            PhysicalTime time_diff = totalTimeSteps*timeStepLength - times.back();

            PhysicalTime time_diff_ratio = time_diff / (entry->first - times.back());
            PhysicalSpeed flux_diff = entry->second - values.back();

            PhysicalSpeed final_flux = values.back() + time_diff_ratio * flux_diff;

            times.push_back(totalTimeSteps*timeStepLength);
            values.push_back(final_flux);
            break;
          }

          times.push_back(entry->first);
          values.push_back(entry->second);
        }
//        densityMin = units->ConvertPressureToLatticeUnits(pMin) / Cs2;
//        densityMax = units->ConvertPressureToLatticeUnits(pMax) / Cs2;

        /* If the time values in the input file end BEFORE the planned end of the simulation, then loop the profile afterwards (using %TimeStepsInInletVelocityProfile). */
        int TimeStepsInStentFluxProfile = times.back() / timeStepLength;

        // Check if last point's value matches the first
        if (values.back() != values.front())
          throw Exception() << "Last point's value does not match the first point's value in "
              << fluxFilePath;

        // extend the table to one past the total time steps, so that the table is valid in the end-state, where the zero indexed time step is equal to the limit.
        fluxTable.resize(totalTimeSteps + 1);
        // Now convert these vectors into arrays using linear interpolation
        for (unsigned int timeStep = 0; timeStep <= totalTimeSteps; timeStep++)
        {
          // The "% TimeStepsInInletVelocityProfile" here is to prevent profile stretching (it will loop instead).
          double point = times.front()
              + (static_cast<double>(timeStep % TimeStepsInStentFluxProfile) / static_cast<double>(totalTimeSteps))
                  * (times.back() - times.front());

          PhysicalSpeed fl = util::NumericalFunctions::LinearInterpolate(times, values, point);

          fluxTable[timeStep] = units->ConvertSpeedToLatticeUnits(fl);
        }

      }

      LatticeSpeed StentFileFlux::GetFlux(const LatticeTimeStep t) const
      {

        if (!useWeightsFromFile)
        {
          // v(r) = vMax (1 - r**2 / a**2)
          // where r is the distance from the centreline

          // Get the max velocity
          LatticeSpeed max = fluxTable[t];

          // Brackets to ensure that the scalar multiplies are done before vector * scalar.
          return max;
        }
        else
        {
          LatticeSpeed max = fluxTable[t];
          
          return max;   
        }

      }

      void StentFileFlux::Initialise(const util::UnitConverter* unitConverter)
      {
        log::Logger::Log<log::Warning, log::OnePerCore>("Initializing fStent.");
        units = unitConverter;

        useWeightsFromFile = false;
        #ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
        useWeightsFromFile = false;
        #endif
      }

    }
  }
}
