#include "lb/boundaries/iolets/InOutLetFile.h"
#include "util/fileutils.h"
#include "util/utilityStructs.h"
#include "SimConfig.h"
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

            TimeValuePair.push_back(tvPair);
          }

          datafile.close();

          std::sort(TimeValuePair.begin(), TimeValuePair.end());

          std::vector<double> time(0);
          std::vector<double> value(0);

          // Must convert into vectors since LinearInterpolate works on a pair of vectors
          for (unsigned int ii = 0; ii < TimeValuePair.size(); ii++)
          {
            time.push_back(TimeValuePair[ii].key);
            value.push_back(TimeValuePair[ii].value);
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
