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

      BoundaryValues::BoundaryValues(geometry::LatticeData::SiteType IOtype,
                                     geometry::LatticeData* iLatDat,
                                     SimConfig* iSimConfig,
                                     SimulationState* iSimState,
                                     util::UnitConverter* iUnits) :
        net::IteratedAction(), mState(iSimState), mSimConfig(iSimConfig), mUnits(iUnits)
      {
        nTotIOlets = (IOtype == geometry::LatticeData::INLET_TYPE
          ? (int) iSimConfig->Inlets.size()
          : (int) iSimConfig->Outlets.size());

        ReadParameters(IOtype);

        mComms = new BoundaryComms(mState, (int) mSimConfig->Inlets.size(), IOtype, iLatDat);
        BCproc = mComms->GetBCProcRank();

        InitialiseBoundaryDensities();

        mComms->SendAndReceive(density);
        mComms->WaitAllComms();
      }

      BoundaryValues::~BoundaryValues()
      {
        delete[] density_avg;
        delete[] density_amp;
        delete[] density_phs;
        delete[] density_min;
        delete[] density_max;
        delete[] filename;

        delete mComms;
      }

      inline bool BoundaryValues::IsCurrentProcTheBCProc()
      {
        return topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc();
      }

      void BoundaryValues::RequestComms()
      {
        if (IsCurrentProcTheBCProc())
        {
          unsigned long time_step = mState->GetTimeStep() % mState->GetTimeStepsPerCycle();
          density = &density_cycle[time_step * nTotIOlets];
        }

        mComms->SendAndReceive(density);
      }

      void BoundaryValues::EndIteration()
      {
        for (int i = 0; i < nTotIOlets; i++)
          CommFinished[i] = false;

        mComms->FinishSend();
      }

      void BoundaryValues::Reset()
      {
        mComms->WaitAllComms();
      }

      void BoundaryValues::ResetPrePeriodChange()
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

      void BoundaryValues::ResetPostPeriodChange()
      {
        // This should occur with the new value of time steps per cycle so need to specify
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
          delete[] density_cycle;
          density_cycle = new distribn_t[hemelb::util::NumericalFunctions::max<int>(1, nTotIOlets)
              * mState->GetTimeStepsPerCycle()];
          density = density_cycle;
        }

        InitialiseBoundaryDensities();

        mComms->SendAndReceive(density);
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

          CommFinished[n] = false;
        }

      }

      void BoundaryValues::allocate()
      {
        if (IsCurrentProcTheBCProc())
        {
          density_cycle = new distribn_t[util::NumericalFunctions::max<int>(1, nTotIOlets)
              * mState->GetTimeStepsPerCycle()];
          density = density_cycle;
        }

        density = new distribn_t[nTotIOlets];
        density_avg = new distribn_t[nTotIOlets];
        density_amp = new distribn_t[nTotIOlets];
        density_phs = new distribn_t[nTotIOlets];
        density_min = new distribn_t[nTotIOlets];
        density_max = new distribn_t[nTotIOlets];
        filename = new std::string[nTotIOlets];
        read_from_file = new bool[nTotIOlets];
        CommFinished = new bool[nTotIOlets];
      }

      void BoundaryValues::FindDensityExtrema()
      {

        // Set initial value in cycle as initial min and mac
        for (int i = 0; i < nTotIOlets; i++)
        {
          density_min[i] = density_cycle[i];
          density_max[i] = density_cycle[i];
        }

        unsigned long offset = 0;

        for (unsigned long i = 0; i < mState->GetTimeStepsPerCycle(); i++)
        {
          for (int j = 0; j < nTotIOlets; j++)
          {
            density_min[j] = util::NumericalFunctions::min(density_min[j],
                                                           density_cycle[offset + j]);
            density_max[j] = util::NumericalFunctions::max(density_max[j],
                                                           density_cycle[offset + j]);
          }

          offset += nTotIOlets;
        }
      }

      distribn_t BoundaryValues::GetBoundaryDensity(const int index)
      {
        if (!CommFinished[index])
        {
          mComms->Wait(index);
          CommFinished[index] = true;
        }

        return density[index];
      }

      distribn_t BoundaryValues::GetDensityMin(int iBoundaryId)
      {
        return density_min[iBoundaryId];
      }

      distribn_t BoundaryValues::GetDensityMax(int iBoundaryId)
      {
        return density_max[iBoundaryId];
      }

    }
  }
}
