#include "lb/boundaries/BoundaryValues.h"
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

        FindBCProcRank();

        std::vector<int> *procsList = new std::vector<int>[nTotIOlets];
        mComms = new BoundaryComms*[nTotIOlets];

        nIOlets = 0;
        iolets = std::vector<int>(0);

        for (int i = 0; i < nTotIOlets; i++)
        {
          bool IOletOnThisProc = IsIOletOnThisProc(IOtype, iLatDat, i);

          procsList[i] = GatherProcList(IOletOnThisProc);

          if (IOletOnThisProc || IsCurrentProcTheBCProc())
          {
            nIOlets++;
            iolets.push_back(i);

            mComms[i] = new BoundaryComms(mState, procsList[i], IOletOnThisProc, BCproc);
          }
        }

        InitialiseBoundaryDensities();

        for (int i = 0; i < nIOlets; i++)
        {
          mComms[iolets[i]]->SendAndReceive(&density[iolets[i]]);
          mComms[iolets[i]]->WaitAllComms();
        }

        // Clear up
        delete[] procsList;
      }

      BoundaryValues::~BoundaryValues()
      {
        if (IsCurrentProcTheBCProc())
        {
          delete[] density_cycle;
          delete[] density_period;
        }
        else
        {
          delete[] density;
        }

        delete[] density_avg;
        delete[] density_amp;
        delete[] density_phs;
        delete[] density_min;
        delete[] density_max;
        delete[] filename;

        for (int i = 0; i < nIOlets; i++)
          delete mComms[iolets[i]];
        delete[] mComms;
      }

      void BoundaryValues::FindBCProcRank()
      {
        proc_t BCrank = 0;

        if (IsCurrentProcTheBCProc())
          BCrank = topology::NetworkTopology::Instance()->GetLocalRank();

        // Since only one proc will update BCrank, the sum of all BCrank is the BCproc
        MPI_Allreduce(&BCrank, &BCproc, 1, hemelb::MpiDataType(BCrank), MPI_SUM, MPI_COMM_WORLD);
      }

      bool BoundaryValues::IsIOletOnThisProc(geometry::LatticeData::SiteType IOtype,
                                              geometry::LatticeData* iLatDat,
                                              int iBoundaryId)
      {
        for (site_t i = 0; i < iLatDat->GetLocalFluidSiteCount(); i++)
          if (iLatDat->GetSiteType(i) == IOtype && iLatDat->GetBoundaryId(i) == iBoundaryId)
            return true;

        return false;
      }

      std::vector<int> BoundaryValues::GatherProcList(bool hasBoundary)
      {
        std::vector<int> procsList(0);

        // This is where the info about whether a proc contains the given inlet/outlet is sent
        // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
        int IOletOnThisProc = hasBoundary; // true if inlet i is on this proc

        // These should be bool, but MPI only supports MPI_INT
        // For each inlet/outlet there is an array of length equal to total number of procs.
        // Each stores true/false value. True if proc of rank equal to the index contains
        // the given inlet/outlet.
        int nTotProcs = topology::NetworkTopology::Instance()->GetProcessorCount();
        int *boolList = new int[nTotProcs];

        MPI_Gather(&IOletOnThisProc,
                   1,
                   hemelb::MpiDataType(IOletOnThisProc),
                   boolList,
                   1,
                   hemelb::MpiDataType(boolList[0]),
                   BCproc,
                   MPI_COMM_WORLD);

        if (IsCurrentProcTheBCProc())
        {
          // Now we have an array for each IOlet with true (1) at indices corresponding to
          // processes that are members of that group. We have to convert this into arrays
          // of ints which store a list of processor ranks.
          for (int j = 0; j < nTotProcs; j++)
          {
            if (boolList[j])
            {
              procsList.push_back(j);
            }
          }
        }

        // Clear up
        delete[] boolList;

        return procsList;
      }

      bool BoundaryValues::IsCurrentProcTheBCProc()
      {
        return topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc();
      }

      void BoundaryValues::RequestComms()
      {
        if (IsCurrentProcTheBCProc())
        {
          UpdateBoundaryDensities();

          unsigned long time_step = (mState->GetTimeStep() - 1) % mState->GetTimeStepsPerCycle();
          density = &density_cycle[time_step * nTotIOlets];
        }

        for (int i = 0; i < nIOlets; i++)
          mComms[iolets[i]]->SendAndReceive(&density[iolets[i]]);
      }

      void BoundaryValues::EndIteration()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          CommFinished[iolets[i]] = false;
          mComms[iolets[i]]->FinishSend();
        }
      }

      void BoundaryValues::ResetPrePeriodChange()
      {
        for (int i = 0; i < nTotIOlets; i++)
        {
          if (!read_from_file[i])
          {
            density_avg[i] = mUnits->ConvertPressureToPhysicalUnits(density_avg[i] * Cs2);
            density_amp[i] = mUnits->ConvertPressureGradToPhysicalUnits(density_amp[i] * Cs2);
            density_min[i] = mUnits->ConvertPressureToPhysicalUnits(density_min[i] * Cs2);
            density_max[i] = mUnits->ConvertPressureToPhysicalUnits(density_max[i] * Cs2);
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
            density_min[i] = mUnits->ConvertPressureToLatticeUnits(density_min[i]) / Cs2;
            density_max[i] = mUnits->ConvertPressureToLatticeUnits(density_max[i]) / Cs2;
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

        for (int i = 0; i < nIOlets; i++)
          mComms[iolets[i]]->SendAndReceive(&density[iolets[i]]);
      }

      void BoundaryValues::Reset()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          mComms[iolets[i]]->WaitAllComms();
        }
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
              if (density_period[i] == 0)
                InitialiseCosCycle(i);
              else
                UpdateCosCycle(i);
            }
          }
        }
      }

      // Should only be called by BCproc
      void BoundaryValues::UpdateBoundaryDensities()
      {
        for (int i = 0; i < nTotIOlets; i++)
          if (density_period[i] != 0)
            if ( (mState->GetTimeStep() - 1) % density_period[i] == 0)
              UpdateCosCycle(i);
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

      void BoundaryValues::UpdateCosCycle(int i)
      {
        double w = 2.0 * PI / (double) mState->GetTimeStepsPerCycle();

        for (unsigned long time_step = mState->GetTimeStep() - 1; time_step
            < (mState->GetTimeStep() - 1) + density_period[i]; time_step++)
        {
          unsigned long tStep = time_step % mState->GetTimeStepsPerCycle();

          density_cycle[tStep * nTotIOlets + i] = density_avg[i] + density_amp[i] * cos(w
              * (double) tStep + density_phs[i]);
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
      void BoundaryValues::SortValuesFromFile(std::vector<double> &time,
                                               std::vector<double> &value)
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

          if (IsCurrentProcTheBCProc())
          {
            density_period[n] = 0;
          }

          filename[n] = lIOlet->PFilePath;
          read_from_file[n] = ! (filename[n] == "");

          if (!read_from_file[n])
          {
            density_avg[n] = mUnits->ConvertPressureToLatticeUnits(lIOlet->PMean) / Cs2;
            density_amp[n] = mUnits->ConvertPressureGradToLatticeUnits(lIOlet->PAmp) / Cs2;
            density_phs[n] = lIOlet->PPhase * DEG_TO_RAD;
          }
          density_min[n] = mUnits->ConvertPressureToLatticeUnits(lIOlet->PMin) / Cs2;
          density_max[n] = mUnits->ConvertPressureToLatticeUnits(lIOlet->PMax) / Cs2;

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

          density_period = new unsigned long[nTotIOlets];
        }
        else
        {
          density = new distribn_t[nTotIOlets];
        }

        density_avg = new distribn_t[nTotIOlets];
        density_amp = new distribn_t[nTotIOlets];
        density_phs = new distribn_t[nTotIOlets];
        density_min = new distribn_t[nTotIOlets];
        density_max = new distribn_t[nTotIOlets];
        filename = new std::string[nTotIOlets];
        read_from_file = new bool[nTotIOlets];
        CommFinished = new bool[nTotIOlets];
      }

      distribn_t BoundaryValues::GetBoundaryDensity(const int index)
      {
        if (!CommFinished[index])
        {
          mComms[index]->Wait();
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
