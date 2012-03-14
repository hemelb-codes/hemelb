#include <algorithm>
#include <fstream>

#include "lb/boundaries/iolets/InOutLetFile.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "util/utilityStructs.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {
        InOutLetFile::InOutLetFile() :
          InOutLet(),densityTable(0)
        {

        }

        void InOutLetFile::DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig)
        {
          iSimConfig->DoIOForFileInOutlet(iParent, iIsLoading, this);
        }

        InOutLet* InOutLetFile::Clone()
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
        void InOutLetFile::CalculateTable(unsigned long total_time_steps)
        {
          // First read in values from file into vectors
          std::vector<util::key_value_pair<double, double> > TimeValuePair(0);

          double timeTemp, valueTemp;

          util::check_file(PressureFilePath.c_str());
          std::ifstream datafile(PressureFilePath.c_str());
          log::Logger::Log<log::Debug, log::OnePerCore>("Reading iolet values from file:");
          while (datafile.good())
          {
            datafile >> timeTemp >> valueTemp;
            log::Logger::Log<log::Debug, log::OnePerCore>("Time: %f Value: %f", timeTemp, valueTemp);
            util::key_value_pair<double, double> tvPair;
            tvPair.key = timeTemp;
            tvPair.value = valueTemp;

            // Don't enter repeat values
            if (TimeValuePair.empty())
            {
              TimeValuePair.push_back(tvPair);
            }
            else
            {
              if (tvPair.key != TimeValuePair[TimeValuePair.size() - 1].key || tvPair.value
                  != TimeValuePair[TimeValuePair.size() - 1].value)
              {
                TimeValuePair.push_back(tvPair);
              }
            }
          }

          datafile.close();

          std::sort(TimeValuePair.begin(), TimeValuePair.end());

          std::vector<double> time(0);
          std::vector<double> value(0);

          // Must convert into vectors since LinearInterpolate works on a pair of vectors
          // Determine min and max pressure on the way
          PressureMinPhysical = TimeValuePair[0].value;
          PressureMaxPhysical = TimeValuePair[0].value;
          for (unsigned int ii = 0; ii < TimeValuePair.size(); ii++)
          {
            PressureMinPhysical = util::NumericalFunctions::min(PressureMinPhysical, TimeValuePair[ii].value);
            PressureMaxPhysical = util::NumericalFunctions::max(PressureMaxPhysical, TimeValuePair[ii].value);
            time.push_back(TimeValuePair[ii].key);
            value.push_back(TimeValuePair[ii].value);
          }

          // Check if last point's value matches the first
          if (value[value.size() - 1] != value[0])
          {
            log::Logger::Log<log::Info, log::OnePerCore>("Last point's value does not match the first point's value in %s\nExiting.",
                                                         PressureFilePath.c_str());
            exit(0);
          }
          // extend the table to one past the total time steps, so that the table is valid in the end-state, where the zero indexed time step is equal to the limit.
          densityTable.resize(total_time_steps+1);
          // Now convert these vectors into arrays using linear interpolation
          for (unsigned int time_step = 0; time_step <= total_time_steps; time_step++)
          {
            double point = time[0] + (static_cast<double>(time_step) / static_cast<double>(total_time_steps)) * (time[time.size() - 1]
                - time[0]);

            double pressure = util::NumericalFunctions::LinearInterpolate(time, value, point);

            densityTable[time_step] = mUnits->ConvertPressureToLatticeUnits(pressure) / Cs2;
          }
        }

      }
    }
  }
}
