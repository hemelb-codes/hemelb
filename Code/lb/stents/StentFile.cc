
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <algorithm>
#include <fstream>

#include "lb/stents/StentFile.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "util/utilityFunctions.h"
#include "util/utilityStructs.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {
      StentFile::StentFile() :
        Stent(), densityTable(0), units(NULL)
      {

      }

      Stent* StentFile::Clone() const
      {
        StentFile* copy = new StentFile(*this);

        return copy;
      }

      StentFile::~StentFile()
      {

      }
      void StentFile::Initialise(const util::UnitConverter* unitConverter)
      {
        units = unitConverter;
      }
      // This reads in a file and interpolates between points to generate a cycle
      // IMPORTANT: to allow reading in data taken at irregular intervals the user
      // needs to make sure that the last point in the file coincides with the first
      // point of a new cycle for a continuous trace.
      void StentFile::CalculateTable(LatticeTimeStep totalTimeSteps, PhysicalTime timeStepLength)
      {
        // First read in values from file
        // Used to be complex code here to keep a vector unique, but this is just achieved by using a map.
        std::map < PhysicalTime, PhysicalDensity > timeValuePairs;

        double timeTemp, valueTemp;

        util::check_file(densityFilePath.c_str());
        std::ifstream datafile(densityFilePath.c_str());
        log::Logger::Log<log::Debug, log::OnePerCore>("Reading stent values from file:");
        while (datafile.good())
        {
          datafile >> timeTemp >> valueTemp;
          log::Logger::Log<log::Trace, log::OnePerCore>("Time: %f Value: %f", timeTemp, valueTemp);
          timeValuePairs[timeTemp] = valueTemp;
        }

        datafile.close();
        // the default iterator for maps traverses in key order, so no sort is needed.

        std::vector<double> times(0);
        std::vector<double> values(0);

        // Must convert into vectors since LinearInterpolate works on a pair of vectors
        // Determine min and max pressure on the way
        PhysicalDensity cMin = timeValuePairs.begin()->second;
        PhysicalDensity cMax = timeValuePairs.begin()->second;
        for (std::map<PhysicalTime, PhysicalDensity>::iterator entry = timeValuePairs.begin(); entry
            != timeValuePairs.end(); entry++)
        {
          /* If the time value stretches beyond the end of the simulation, then insert an interpolated end value and exit the loop. */
          if(entry->first > totalTimeSteps*timeStepLength) {

            PhysicalTime time_diff = totalTimeSteps*timeStepLength - times.back();

            PhysicalTime time_diff_ratio = time_diff / (entry->first - times.back());
            PhysicalDensity c_diff = entry->second - values.back();

            PhysicalDensity final_concentration = values.back() + time_diff_ratio * c_diff;

            times.push_back(totalTimeSteps*timeStepLength);
            cMin = util::NumericalFunctions::min(cMin, final_concentration);
            cMax = util::NumericalFunctions::max(cMax, final_concentration);
            values.push_back(final_concentration);
            break;
          }

          cMin = util::NumericalFunctions::min(cMin, entry->second);
          cMax = util::NumericalFunctions::max(cMax, entry->second);
          times.push_back(entry->first);
          values.push_back(entry->second);
        }
        densityMin = units->ConvertDensityToLatticeUnits(cMin);
        densityMax = units->ConvertDensityToLatticeUnits(cMax);

        /* If the time values in the input file end BEFORE the planned end of the simulation, then loop the profile afterwards (using %TimeStepsInInletPressureProfile). */
        int TimeStepsInStentConcentrationProfile = times.back() / timeStepLength;
        //throw Exception() << "Finding Time steps: " << times.back() << " " << timeStepLength << " " << TimeStepsInInletPressureProfile;

        // extend the table to one past the total time steps, so that the table is valid in the end-state, where the zero indexed time step is equal to the limit.
        densityTable.resize(totalTimeSteps + 1);
        // Now convert these vectors into arrays using linear interpolation
        for (unsigned int timeStep = 0; timeStep <= totalTimeSteps; timeStep++)
        {
          double point = times.front() + (static_cast<double> (timeStep)
              / static_cast<double> (totalTimeSteps)) * (times.back() - times.front());

          double concentration = util::NumericalFunctions::LinearInterpolate(times, values, point);

          densityTable[timeStep] = units->ConvertDensityToLatticeUnits(concentration);
        }
      }

    }
  }
}
