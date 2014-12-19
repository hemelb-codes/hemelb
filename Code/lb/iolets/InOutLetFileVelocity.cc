//
// Copyright (C) University College London, 2007-2013, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "lb/iolets/InOutLetFileVelocity.h"
#include <algorithm>
#include <fstream>
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
      InOutLetFileVelocity::InOutLetFileVelocity() :
          units(NULL)
      {
      }

      InOutLet* InOutLetFileVelocity::Clone() const
      {
        InOutLet* copy = new InOutLetFileVelocity(*this);
        return copy;
      }

      void InOutLetFileVelocity::CalculateTable(LatticeTimeStep totalTimeSteps)
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
//          pMin = util::NumericalFunctions::min(pMin, entry->second);
//          pMax = util::NumericalFunctions::max(pMax, entry->second);
          times.push_back(entry->first);
          values.push_back(entry->second);
        }
//        densityMin = units->ConvertPressureToLatticeUnits(pMin) / Cs2;
//        densityMax = units->ConvertPressureToLatticeUnits(pMax) / Cs2;

        // Check if last point's value matches the first
        if (values.back() != values.front())
          throw Exception() << "Last point's value does not match the first point's value in "
              << velocityFilePath;

        // extend the table to one past the total time steps, so that the table is valid in the end-state, where the zero indexed time step is equal to the limit.
        velocityTable.resize(totalTimeSteps + 1);
        // Now convert these vectors into arrays using linear interpolation
        for (unsigned int timeStep = 0; timeStep <= totalTimeSteps; timeStep++)
        {
          double point = times.front()
              + (static_cast<double>(timeStep) / static_cast<double>(totalTimeSteps))
                  * (times.back() - times.front());

          PhysicalSpeed vel = util::NumericalFunctions::LinearInterpolate(times, values, point);

          velocityTable[timeStep] = units->ConvertVelocityToLatticeUnits(vel);
        }

      }

      LatticeVelocity InOutLetFileVelocity::GetVelocity(const LatticePosition& x,
                                                        const LatticeTimeStep t) const
      {

        if (!useWeightsFromFile)
        {
          // v(r) = vMax (1 - r**2 / a**2)
          // where r is the distance from the centreline
          LatticePosition displ = x - position;
          LatticeDistance z = displ.Dot(normal);
          Dimensionless rSqOverASq = (displ.GetMagnitudeSquared() - z * z) / (radius * radius);
          assert(rSqOverASq <= 1.0);

          // Get the max velocity
          LatticeSpeed max = velocityTable[t];

          // Brackets to ensure that the scalar multiplies are done before vector * scalar.
          return normal * (max * (1. - rSqOverASq));
        }
        else
        {
          std::vector<int> xyz;
          xyz.push_back((int) (x.x + 0.5));
          xyz.push_back((int) (x.y + 0.5));
          xyz.push_back((int) (x.z + 0.5));

          int half[3];
          int valid_count = 1;
          if (x.x - int(x.x) > 0.1)
          {
            half[0] = 1;
            valid_count *= 2;
          }
          if (x.y - int(x.y) > 0.1)
          {
            half[1] = 1;
            valid_count *= 2;
          }
          if (x.z - int(x.z) > 0.1)
          {
            half[2] = 1;
            valid_count *= 2;
          }

          double v[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
          double v_tot = 0.0;
          for (int i = 0; i < half[0]; i++)
          {
            for (int j = 0; i < half[1]; i++)
            {
              for (int k = 0; i < half[2]; i++)
              {
                try
                {
                  v[i * 4 + j * 2 + k] = weights_table.at(xyz);
                  v_tot += v[i * 4 + j * 2 + k];
                }
                catch (int e)
                {
                  v[i * 4 + j * 2 + k] = 0;
                  valid_count--; //no matching entry, so we'll drop the valid_count.
                }
              }
            }
          }
          return v_tot / valid_count;
        }

      }

      void InOutLetFileVelocity::Initialise(const util::UnitConverter* unitConverter)
      {
        log::Logger::Log<log::Warning, log::OnePerCore>("Initializing vInlet.");
        units = unitConverter;


        useWeightsFromFile = true;

        //if the new velocity approximation is enabled, then we want to create a lookup table here.
        const std::string in_name = velocityFilePath + ".weights.txt";

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
        while (std::getline(myfile, input_line))
        {
          int x, y, z;
          double v;
          myfile >> x >> y >> z >> v;

          std::vector<int> xyz;
          xyz.push_back(x);
          xyz.push_back(y);
          xyz.push_back(z);
          weights_table[xyz] = v;

          log::Logger::Log<log::Warning, log::OnePerCore>("%lld %lld %lld %f",
                                                          x,
                                                          y,
                                                          z,
                                                          weights_table[xyz]);
        }
      }

    }
  }
}
