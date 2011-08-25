#include "lb/boundaries/InOutLetFile.h"
#include "util/fileutils.h"
#include <fstream>
#include <algorithm>

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      // This is used to sort read in values
      // Allows sorting time-value pairs using the standard library sort
      struct time_value_pair
      {
        public:
          double time;
          double value;

          bool operator<(const time_value_pair other_time_value_pair) const
          {
            return time < other_time_value_pair.time;
          }
      };

      InOutLetFile::InOutLetFile(unsigned long iPeriod,
                                 double iPMin,
                                 double iPMax,
                                 util::Vector3D iPosition,
                                 util::Vector3D iNormal,
                                 util::UnitConverter* iUnits,
                                 std::string iPFilePath) :
        InOutLet(iPeriod, iPMin, iPMax, iPosition, iNormal, iUnits), PressureFilePath(iPFilePath)
      {

      }

      InOutLetFile::~InOutLetFile()
      {

      }

      void InOutLetFile::CalculateCycle(std::vector<distribn_t> &density_cycle,
                                        const SimulationState *iState)
      {
        // First read in values from file into vectors
        std::vector<time_value_pair> TimeValuePair(0);

        double timeTemp, valueTemp;

        util::check_file(PressureFilePath.c_str());
        std::ifstream datafile(PressureFilePath.c_str());

        while (datafile.good())
        {
          datafile >> timeTemp >> valueTemp;
          time_value_pair tvPair;
          tvPair.time = timeTemp;
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
          time.push_back(TimeValuePair[ii].time);
          value.push_back(TimeValuePair[ii].value);
        }

        // Now convert these vectors into arrays using linear interpolation
        for (unsigned int time_step = 0; time_step < density_cycle.size(); time_step++)
        {
          double point = time[0] + (double) time_step / (double) density_cycle.size()
              * (time[time.size() - 1] - time[0]);

          double pressure = util::NumericalFunctions::LinearInterpolate(time, value, point);

          density_cycle[time_step] = mUnits->ConvertPressureToLatticeUnits(pressure) / Cs2;
        }
      }

    }
  }
}
