#include "lb/boundaries/BoundaryValues.h"
#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include "util/fileutils.h"
#include <fstream>
#include <math.h>

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      BoundaryValues::BoundaryValues(BoundaryComms* iComms,
                                     geometry::LatticeData::SiteType IOtype,
                                     geometry::LatticeData* iLatDat,
                                     SimConfig* iSimConfig,
                                     SimulationState* iSimState,
                                     util::UnitConverter* iUnits) :
        mState(iSimState), mSimConfig(iSimConfig), mUnits(iUnits)
      {
        proc_t BCrank = 0;

        if (IsCurrentProcTheBCProc())
          BCrank = topology::NetworkTopology::Instance()->GetLocalRank();

        // Since only one proc will update BCrank, the sum of all BCrank is the BCproc
        MPI_Allreduce(&BCrank, &BCproc, 1, hemelb::MpiDataType(BCrank), MPI_SUM, MPI_COMM_WORLD);

        nTotIOlets = (IOtype == geometry::LatticeData::INLET_TYPE
          ? (int) iSimConfig->Inlets.size()
          : (int) iSimConfig->Outlets.size());

        density_cycle = std::vector<distribn_t>(0);

        ReadParameters(IOtype);

        InitialiseBoundaryDensities();

        iComms->Initialise(IOtype, iLatDat, &density_cycle);
      }

      BoundaryValues::~BoundaryValues()
      {
        delete[] density_avg;
        delete[] density_amp;
        delete[] density_phs;
        delete[] density_min;
        delete[] density_max;
        delete[] filename;
      }

      inline bool BoundaryValues::IsCurrentProcTheBCProc()
      {
        return topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc();
      }

      void BoundaryValues::InitialiseBoundaryDensities()
      {
        if (IsCurrentProcTheBCProc())
        {
          for (int i = 0; i < nTotIOlets; i++)
          {
            if (read_from_file[i])
            {
              InitialiseFromFile(i);
            }
            else
            {
              InitialiseCosCycle(i);
            }
          }

          FindDensityExtrema();
        }

        MPI_Bcast(density_min,
                  nTotIOlets,
                  hemelb::MpiDataType(density_min[0]),
                  BCproc,
                  MPI_COMM_WORLD);

        MPI_Bcast(density_max,
                  nTotIOlets,
                  hemelb::MpiDataType(density_max[0]),
                  BCproc,
                  MPI_COMM_WORLD);
      }

      void BoundaryValues::InitialiseCosCycle(int i)
      {
        double w = 2.0 * PI / (double) mState->GetTimeStepsPerCycle();

        for (unsigned long time_step = 0; time_step < mState->GetTimeStepsPerCycle(); time_step++)
        {
          density_cycle[time_step * nTotIOlets + i] = density_avg[i] + density_amp[i] * cos(w
              * (double) time_step + density_phs[i]);
        }
      }

      void BoundaryValues::InitialiseFromFile(int i)
      {
        // First read in values from file into vectors

        std::vector<double> time(0);
        std::vector<double> value(0);

        double timeTemp, valueTemp;

        util::check_file(filename[i].c_str());
        std::ifstream datafile(filename[i].c_str());

        while (datafile.good())
        {
          datafile >> timeTemp >> valueTemp;
          time.push_back(timeTemp);
          value.push_back(valueTemp);
        }

        datafile.close();

        SortValuesFromFile(time, value);

        // Now convert these vectors into arrays using linear interpolation

        for (unsigned long time_step = 0; time_step < mState->GetTimeStepsPerCycle(); time_step++)
        {
          double point = time[0] + (double) time_step / (double) mState->GetTimeStepsPerCycle()
              * (time[time.size() - 1] - time[0]);

          double pressure = util::NumericalFunctions::LinearInterpolate(time, value, point);

          density_cycle[time_step * nTotIOlets + i]
              = mUnits->ConvertPressureToLatticeUnits(pressure) / Cs2;
        }
      }

      // Sorts both the time and value vector so that time values are in incrementing order
      // Uses bubble sort as it is used only when initialising and reseting
      void BoundaryValues::SortValuesFromFile(std::vector<double> &time, std::vector<double> &value)
      {
        bool swapped = true;

        while (swapped)
        {
          swapped = false;
          for (int i = 0; i < (int) time.size() - 1; i++)
          {
            if (time[i] > time[i + 1])
            {
              double temp = time[i];
              time[i] = time[i + 1];
              time[i + 1] = temp;

              temp = value[i];
              value[i] = value[i + 1];
              value[i + 1] = temp;

              swapped = true;
            }
          }
        }
      }

      void BoundaryValues::ReadParameters(geometry::LatticeData::SiteType IOtype)
      {
        allocate();

        for (int n = 0; n < nTotIOlets; n++)
        {
          hemelb::SimConfig::InOutLet *lIOlet = (IOtype == geometry::LatticeData::INLET_TYPE
            ? &mSimConfig->Inlets[n]
            : &mSimConfig->Outlets[n]);

          density_avg[n] = mUnits->ConvertPressureToLatticeUnits(lIOlet->PMean) / Cs2;
          density_amp[n] = mUnits->ConvertPressureGradToLatticeUnits(lIOlet->PAmp) / Cs2;
          density_phs[n] = lIOlet->PPhase * DEG_TO_RAD;
          filename[n] = lIOlet->PFilePath;

          if (filename[n] == "")
            read_from_file[n] = false;
          else
            read_from_file[n] = true;
        }

      }

      void BoundaryValues::allocate()
      {
        if (IsCurrentProcTheBCProc())
        {
          density_cycle.resize(util::NumericalFunctions::max<int>(1, nTotIOlets)
              * mState->GetTimeStepsPerCycle());
        }

        density_avg = new distribn_t[nTotIOlets];
        density_amp = new distribn_t[nTotIOlets];
        density_phs = new distribn_t[nTotIOlets];
        density_min = new distribn_t[nTotIOlets];
        density_max = new distribn_t[nTotIOlets];
        filename = new std::string[nTotIOlets];
        read_from_file = new bool[nTotIOlets];
      }

      void BoundaryValues::FindDensityExtrema()
      {
        for (int i = 0; i < nTotIOlets; i++)
        {
          density_min[i] = density_cycle[0];
          density_max[i] = density_cycle[0];

          for (unsigned int j = 0; j < density_cycle.size(); j++)
          {
            density_min[i] = util::NumericalFunctions::min(density_min[i], density_cycle[j]);
            density_max[i] = util::NumericalFunctions::max(density_max[i], density_cycle[j]);
          }
        }
      }

      distribn_t BoundaryValues::GetDensityMin(int iBoundaryId)
      {
        return density_min[iBoundaryId];
      }

      distribn_t BoundaryValues::GetDensityMax(int iBoundaryId)
      {
        return density_max[iBoundaryId];
      }

      void BoundaryValues::ResetPrePeriodDoubling()
      {
        for (int i = 0; i < nTotIOlets; i++)
        {
          if (!read_from_file[i])
          {
            density_avg[i] = mUnits->ConvertPressureToPhysicalUnits(density_avg[i] * Cs2);
            density_amp[i] = mUnits->ConvertPressureGradToPhysicalUnits(density_amp[i] * Cs2);
          }
        }
      }

      void BoundaryValues::ResetPostPeriodDoubling()
      {
        for (int i = 0; i < nTotIOlets; i++)
        {
          if (!read_from_file[i])
          {
            density_avg[i] = mUnits->ConvertPressureToLatticeUnits(density_avg[i]) / Cs2;
            density_amp[i] = mUnits->ConvertPressureGradToLatticeUnits(density_amp[i]) / Cs2;
          }
        }

        if (IsCurrentProcTheBCProc())
        {
          density_cycle.resize(hemelb::util::NumericalFunctions::max<int>(1, nTotIOlets)
              * mState->GetTimeStepsPerCycle());
        }

        InitialiseBoundaryDensities();
      }

    }
  }
}
