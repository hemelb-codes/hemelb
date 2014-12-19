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

      int InOutLetFileVelocity::FindWeightIndex(std::string id) const
      {
        for (unsigned int i = 0; i < weight_ids.size(); i++)
        {
          if ( (weight_ids[i]).compare(id) == 0)
          {
            return i;
          }
        }
        //log::Logger::Log<log::Warning, log::OnePerCore>("Error: could not find proper bin: %s",
        //                                                id.c_str());
        return -1;
      }

      LatticeVelocity InOutLetFileVelocity::GetVelocity(const LatticePosition& x,
                                                        const LatticeTimeStep t) const
      {
        // Get the max velocity
        LatticeSpeed max = velocityTable[t];

        if (!useWeightsFromFile)
        {
          // v(r) = vMax (1 - r**2 / a**2)
          // where r is the distance from the centreline
          LatticePosition displ = x - position;
          LatticeDistance z = displ.Dot(normal);
          Dimensionless rSqOverASq = (displ.GetMagnitudeSquared() - z * z) / (radius * radius);
          assert(rSqOverASq <= 1.0);

          // Brackets to ensure that the scalar multiplies are done before vector * scalar.
          return normal * (max * (1. - rSqOverASq));
        }
        else
        {
          /* Simplistic and error-prone conversion
           * TODO: incorporate correct support for .5 position values
           * by interpolating between entries. */
          int xyz[3];
          int xhalf = 0;
          int yhalf = 0;
          int zhalf = 0;
          if (x.x - int(x.x) > 0.1)
          {
            xhalf = 1;
          }
          if (x.y - int(x.y) > 0.1)
          {
            yhalf = 1;
          }
          if (x.z - int(x.z) > 0.1)
          {
            zhalf = 1;
          }

          xyz[0] = (int) (x.x + 0.5);
          xyz[1] = (int) (x.y + 0.5);
          xyz[2] = (int) (x.z + 0.5);

          log::Logger::Log<log::Warning, log::OnePerCore>("Weight index search: %f %f %f",
                                                          x.x,
                                                          x.y,
                                                          x.z);

          for (int i = 0; i < xhalf; i++)
          {
            for (int j = 0; i < yhalf; i++)
            {
              for (int k = 0; i < zhalf; i++)
              {
                std::ostringstream xyz_str;
                xyz_str << xyz[0]+i << " " << xyz[1]+j << " " << xyz[2]+k;
                /* This is very slow, and we can speed this up by building a good/better hash table-like structure. */
                int index = this->FindWeightIndex(xyz_str.str());
              }
            }
          }

          log::Logger::Log<log::Warning, log::OnePerCore>("Weight index found: %d %d",
                                                          index,
                                                          weight_ids.size());

          LatticeSpeed vw = weights[index];
          return normal * max * vw;
        }
      }

      //void PopulateWeightsTableBasedOnCircularInflow() {
      //  weights_table.insert(std::pair<std::vector<int>, PhysicalVelocity>());
      //}

      void InOutLetFileVelocity::Initialise(const util::UnitConverter* unitConverter)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Initializing vInlet.");
        units = unitConverter;

        //if the new velocity approximation is enabled, then we want to create a lookup table here.
        const std::string in_name = velocityFilePath + ".weights.txt";

        /* Load and read file. */
        std::fstream myfile;
        myfile.open(in_name.c_str(), std::ios_base::in);

        std::string input_line;
        /* input files are in ASCII, in format:
         *
         * coord_x coord_y coord_z weights_value
         *
         * */

        /* Temporary for testing. */
        useWeightsFromFile = true;

        log::Logger::Log<log::Warning, log::OnePerCore>("Reading weights file: %s",
                                                        in_name.c_str());
        int x, y, z;
        double v;
        while (myfile >> x >> y >> z >> v)
        {
          std::ostringstream xyz;
          xyz << x << " " << y << " " << z;
          weight_ids.push_back(xyz.str());
          weights.push_back(v);

          log::Logger::Log<log::Debug, log::OnePerCore>("Read %d %d %d %f",
                                                        x,
                                                        y,
                                                        z,
                                                        weights[weights.size() - 1]);
        }
      }

    }
  }
}
