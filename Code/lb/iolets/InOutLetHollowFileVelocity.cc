
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLetHollowFileVelocity.h"
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
    namespace iolets
    {
      InOutLetHollowFileVelocity::InOutLetHollowFileVelocity() :
          units(NULL)
      {
      }

      InOutLet* InOutLetHollowFileVelocity::Clone() const
      {
        InOutLet* copy = new InOutLetHollowFileVelocity(*this);
        return copy;
      }

      void InOutLetHollowFileVelocity::CalculateTable(LatticeTimeStep totalTimeSteps, PhysicalTime timeStepLength)
      {
        // First read in values from file
        // Used to be complex code here to keep a vector unique, but this is just achieved by using a map.
        std::map<PhysicalTime, PhysicalSpeed> timeValuePairs;

        double timeTemp, valueTemp;

        util::check_file(velocityFilePath.c_str());
        std::ifstream datafile(velocityFilePath.c_str());
        log::Logger::Log<log::Debug, log::OnePerCore>("Reading iolet values from file:");
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
            PhysicalSpeed vel_diff = entry->second - values.back();

            PhysicalSpeed final_velocity = values.back() + time_diff_ratio * vel_diff;

            times.push_back(totalTimeSteps*timeStepLength);
            values.push_back(final_velocity);
            break;
          }

          times.push_back(entry->first);
          values.push_back(entry->second);
        }
//        densityMin = units->ConvertPressureToLatticeUnits(pMin) / Cs2;
//        densityMax = units->ConvertPressureToLatticeUnits(pMax) / Cs2;

        /* If the time values in the input file end BEFORE the planned end of the simulation, then loop the profile afterwards (using %TimeStepsInInletVelocityProfile). */
        int TimeStepsInInletVelocityProfile = times.back() / timeStepLength;

        // Check if last point's value matches the first
        //if (values.back() != values.front())
          //throw Exception() << "Last point's value does not match the first point's value in "
              //<< velocityFilePath;

        // extend the table to one past the total time steps, so that the table is valid in the end-state, where the zero indexed time step is equal to the limit.
        velocityTable.resize(totalTimeSteps + 1);
        // Now convert these vectors into arrays using linear interpolation
        for (unsigned int timeStep = 0; timeStep <= totalTimeSteps; timeStep++)
        {
          // The "% TimeStepsInInletVelocityProfile" here is to prevent profile stretching (it will loop instead).
          double point = times.front()
              + (static_cast<double>(timeStep % TimeStepsInInletVelocityProfile) / static_cast<double>(totalTimeSteps))
                  * (times.back() - times.front());

          PhysicalSpeed vel = util::NumericalFunctions::LinearInterpolate(times, values, point);

          velocityTable[timeStep] = units->ConvertVelocityToLatticeUnits(vel);
        }

      }

      LatticeVelocity InOutLetHollowFileVelocity::GetVelocity(const LatticePosition& x,
                                                        const LatticeTimeStep t) const
      {

        if (!useWeightsFromFile)
        {
          // v(r) = (-4 * vMax / ((r_o - r_i) ^ 2)) * (r - r_o) * (r - r_i)
          // where r is the distance from the centreline
          LatticePosition displ = x - position;
          LatticeDistance z = displ.Dot(normal);
          LatticeDistance r_distance = sqrt(displ.GetMagnitudeSquared() - z * z);
          LatticeDistance r_o = radius;
          LatticeDistance r_i = innerRadius;
          
          // Get the max velocity
          LatticeSpeed max = velocityTable[t];

          // Brackets to ensure that the scalar multiplies are done before vector * scalar.
          return normal * (((-4.0 * max) / ((r_o - r_i) * (r_o - r_i))) * (r_distance - r_o) * (r_distance - r_i));
        }
        else
        {
          log::Logger::Log<log::Error, log::Singleton>("Velocity weights file cannot be used for hollow file velocity iolet boundary condition.");
        }

      }

      void InOutLetHollowFileVelocity::Initialise(const util::UnitConverter* unitConverter)
      {
        log::Logger::Log<log::Warning, log::OnePerCore>("Initializing vInlet.");
        units = unitConverter;

        useWeightsFromFile = false;
        #ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
        useWeightsFromFile = true;
        #endif

        if(useWeightsFromFile) {
          //if the new velocity approximation is enabled, then we want to create a lookup table here.
          const std::string in_name = velocityFilePath + ".weights.txt";
          util::check_file(in_name.c_str());

          /* Load and read file. */
          std::fstream myfile;
          myfile.open(in_name.c_str(), std::ios_base::in);
          log::Logger::Log<log::Warning, log::OnePerCore>("Loading weights file: %s",
                                                        in_name.c_str());

          std::string input_line;
          /* input files are in ASCII, in format:
           *
           * coord_x coord_y coord_z weights_value
           *
           * */
          while (myfile.good()) //(std::getline(myfile, input_line))
          {
            int x, y, z;
            double v;
            myfile >> x >> y >> z >> v;

            std::vector<int> xyz;
            xyz.push_back(x);
            xyz.push_back(y);
            xyz.push_back(z);
            weights_table[xyz] = v;

            log::Logger::Log<log::Trace, log::OnePerCore>("%lld %lld %lld %f",
            x,
            y,
            z,
            weights_table[xyz]);
          }
          myfile.close();
        }
      }

    }
  }
}
