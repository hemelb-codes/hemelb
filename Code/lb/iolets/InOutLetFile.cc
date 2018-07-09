
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <algorithm>
#include <fstream>

#include "lb/iolets/InOutLetFile.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "util/utilityFunctions.h"
#include "util/utilityStructs.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      InOutLetFile::InOutLetFile() :
        InOutLet(), densityTable(0), units(NULL)
      {

      }

      InOutLet* InOutLetFile::Clone() const
      {
        InOutLetFile* copy = new InOutLetFile(*this);

        return copy;
      }

      InOutLetFile::~InOutLetFile()
      {

      }
      void InOutLetFile::Initialise(const util::UnitConverter* unitConverter)
      {
        units = unitConverter;
      }
      // This reads in a file and interpolates between points to generate a cycle
      // IMPORTANT: to allow reading in data taken at irregular intervals the user
      // needs to make sure that the last point in the file coincides with the first
      // point of a new cycle for a continuous trace.
      void InOutLetFile::CalculateTable(LatticeTimeStep totalTimeSteps, PhysicalTime timeStepLength)
      {
        // First read in values from file
        // Used to be complex code here to keep a vector unique, but this is just achieved by using a map.
        std::map < PhysicalTime, PhysicalPressure > timeValuePairs;

        double timeTemp, valueTemp;

        util::check_file(pressureFilePath.c_str());
        std::ifstream datafile(pressureFilePath.c_str());
        log::Logger::Log<log::Debug, log::OnePerCore>("Reading iolet values from file:");
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
        PhysicalPressure pMin = timeValuePairs.begin()->second;
        PhysicalPressure pMax = timeValuePairs.begin()->second;
        for (std::map<PhysicalTime, PhysicalPressure>::iterator entry = timeValuePairs.begin(); entry
            != timeValuePairs.end(); entry++)
        {
          /* If the time value stretches beyond the end of the simulation, then insert an interpolated end value and exit the loop. */
          if(entry->first > totalTimeSteps*timeStepLength) {

            PhysicalTime time_diff = totalTimeSteps*timeStepLength - times.back();

            PhysicalTime time_diff_ratio = time_diff / (entry->first - times.back());
            PhysicalPressure pres_diff = entry->second - values.back();

            PhysicalSpeed final_pressure = values.back() + time_diff_ratio * pres_diff;

            times.push_back(totalTimeSteps*timeStepLength);
            pMin = util::NumericalFunctions::min(pMin, final_pressure);
            pMax = util::NumericalFunctions::max(pMax, final_pressure);
            values.push_back(final_pressure);
            break;
          }

          pMin = util::NumericalFunctions::min(pMin, entry->second);
          pMax = util::NumericalFunctions::max(pMax, entry->second);
          times.push_back(entry->first);
          values.push_back(entry->second);
        }
        densityMin = units->ConvertPressureToLatticeUnits(pMin) / Cs2;
        densityMax = units->ConvertPressureToLatticeUnits(pMax) / Cs2;

        // Check if last point's value matches the first
        if (values.back() != values.front())
          throw Exception() << "Last point's value does not match the first point's value in " <<pressureFilePath;

        /* If the time values in the input file end BEFORE the planned end of the simulation, then loop the profile afterwards (using %TimeStepsInInletPressureProfile). */
        int TimeStepsInInletPressureProfile = times.back() / timeStepLength;
        //throw Exception() << "Finding Time steps: " << times.back() << " " << timeStepLength << " " << TimeStepsInInletPressureProfile;

        // extend the table to one past the total time steps, so that the table is valid in the end-state, where the zero indexed time step is equal to the limit.
        densityTable.resize(totalTimeSteps + 1);
        // Now convert these vectors into arrays using linear interpolation
        for (unsigned int timeStep = 0; timeStep <= totalTimeSteps; timeStep++)
        {
          double point = times.front() + (static_cast<double> (timeStep)
              / static_cast<double> (totalTimeSteps)) * (times.back() - times.front());

          double pressure = util::NumericalFunctions::LinearInterpolate(times, values, point);

          densityTable[timeStep] = units->ConvertPressureToLatticeUnits(pressure) / Cs2;
        }
      }

    }
  }
}
