// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
        InOutLet(), densityTable(0)
      {

      }

      void InOutLetFile::DoIO(TiXmlElement *parent, bool isLoading,
                              configuration::SimConfig* simConfig)
      {
        simConfig->DoIOForFileInOutlet(parent, isLoading, this);
      }

      InOutLet* InOutLetFile::Clone() const
      {
        InOutLetFile* copy = new InOutLetFile(*this);

        return copy;
      }

      InOutLetFile::~InOutLetFile()
      {

      }

      // This reads in a file and interpolates between points to generate a cycle
      // IMPORTANT: to allow reading in data taken at irregular intervals the user
      // needs to make sure that the last point in the file coincides with the first
      // point of a new cycle for a continuous trace.
      void InOutLetFile::CalculateTable(LatticeTime totalTimeSteps)
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
        pressureMinPhysical = timeValuePairs.begin()->second;
        pressureMaxPhysical = timeValuePairs.begin()->second;
        for (std::map<PhysicalTime, PhysicalPressure>::iterator entry = timeValuePairs.begin(); entry
            != timeValuePairs.end(); entry++)
        {
          pressureMinPhysical = util::NumericalFunctions::min(pressureMinPhysical, entry->second);
          pressureMaxPhysical = util::NumericalFunctions::max(pressureMaxPhysical, entry->second);
          times.push_back(entry->first);
          values.push_back(entry->second);
        }

        // Check if last point's value matches the first
        if (values.back() != values.front())
        {
          log::Logger::Log<log::Critical, log::OnePerCore>("Last point's value does not match the first point's value in %s\nExiting.",
                                                           pressureFilePath.c_str());
          exit(0);
        }
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
