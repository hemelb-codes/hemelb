#include "lb/boundaries/iolets/InOutLetFile.h"
#include "util/fileutils.h"
#include "util/utilityStructs.h"
#include "SimConfig.h"
#include "log/Logger.h"
#include <fstream>
#include <algorithm>

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {
        InOutLetFile::InOutLetFile() :
          InOutLetCycle<0, false> ()
        {

        }

        void InOutLetFile::DoIO(TiXmlElement *iParent, bool iIsLoading, SimConfig* iSimConfig)
        {
          iSimConfig->DoIO(iParent, iIsLoading, this);
        }

        InOutLet* InOutLetFile::Clone()
        {
          InOutLetFile* copy = new InOutLetFile();
          copy->PressureMinPhysical = this->PressureMinPhysical;
          copy->PressureMaxPhysical = this->PressureMaxPhysical;
          copy->Position = this->Position;
          copy->Normal = this->Normal;
          copy->PressureFilePath = this->PressureFilePath;

          return copy;
        }

        InOutLetFile::~InOutLetFile()
        {

        }

        // This reads in a file and interpolates between points to generate a cycle
        // IMPORTANT: to allow reading in data taken at irregular intervals the user
        // needs to make sure that the last point in the file coincides with the first
        // point of a new cycle for a continuous trace.
        void InOutLetFile::CalculateCycle(std::vector<distribn_t> &density_cycle,
                                          const SimulationState *iState)
        {
          // First read in values from file into vectors
          std::vector<util::key_value_pair<double, double> > TimeValuePair(0);

          double timeTemp, valueTemp;

          util::check_file(PressureFilePath.c_str());
          std::ifstream datafile(PressureFilePath.c_str());

          while (datafile.good())
          {
            datafile >> timeTemp >> valueTemp;
            util::key_value_pair<double, double> tvPair;
            tvPair.key = timeTemp;
            tvPair.value = valueTemp;

            // Don't enter repeat values
            if (TimeValuePair.size() == 0)
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

          // Now convert these vectors into arrays using linear interpolation
          for (unsigned int time_step = 0; time_step < density_cycle.size(); time_step++)
          {
            double point = time[0] + ((double) time_step / (double) density_cycle.size())
                * (time[time.size() - 1] - time[0]);

            double pressure = util::NumericalFunctions::LinearInterpolate(time, value, point);

            density_cycle[time_step] = util::UnitConverter::ConvertPressureToLatticeUnits(pressure)
                / Cs2;
          }
        }

      }
    }
  }
}
