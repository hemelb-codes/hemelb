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
          //log::Logger::Log<log::Warning, log::OnePerCore>("%f %f %f", x.x, x.y, x.z);
          std::vector<int> xyz;
          xyz.push_back(0);
          xyz.push_back(0);
          xyz.push_back(0);

          int whole[3] = {1,1,1};
          int xint = int(x.x - 0.9999999);
          int yint = int(x.y - 0.9999999);
          int zint = int(x.z - 0.9999999);

          if (x.x - ((float) xint) > 0.1)
          {
            whole[0] = 0;
            xint += 1;
          }
          if (x.y - ((float) yint) > 0.1)
          {
            whole[1] = 0;
            yint += 1;
          }
          if (x.z - ((float) zint) > 0.1)
          {
            whole[2] = 0;
            zint += 1;
          }

          //log::Logger::Log<log::Warning, log::OnePerCore>("%f %f", x.x, ((float) xint));

          std::vector<double> v;
          v.reserve(27);
          std::vector<bool> dist_flag;
          dist_flag.reserve(27);

          int sizes[3] = {whole[0]+2, whole[1]+2, whole[2]+2};
          int v_tot = 0;

          for (int i = 0; i < sizes[0]; i++)
          {
            for (int j = 0; j < sizes[1]; j++)
            {

              for (int k = 0; k < sizes[2]; k++)
              {

                xyz[0] = xint + i;
                xyz[1] = yint + j;
                xyz[2] = zint + k;

                if(weights_table.count(xyz)>0)
                {
                  int is_dist = (whole[0] * (1-(i%2))) + (whole[1] * (1-(j%2))) + (whole[2] * (1-(k%2)));
                  v.push_back(weights_table.at(xyz));
                  dist_flag.push_back(is_dist);
                  v_tot += weights_table.at(xyz);
                  //log::Logger::Log<log::Warning, log::OnePerCore>("%d %d %d OK.", xyz[0], xyz[1], xyz[2]);
                }
              }
            }
          }
          //log::Logger::Log<log::Warning, log::Singleton>("%f %f %f %d",
          //                                                x.x,
          //                                                x.y,
          //                                                x.z, valid_count);
         
          //if(v.size() == 0) {
          /* local interpolation did not work, so we adopt a wider range. */ 

          //}
 
          return v_tot / v.size();
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

          /*log::Logger::Log<log::Warning, log::OnePerCore>("%lld %lld %lld %f",
           x,
           y,
           z,
           weights_table[xyz]);*/
        }
      }

    }
  }
}
