#include "lb/BoundaryComms.h"
#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include <math.h>

namespace hemelb
{
  namespace lb
  {

    BoundaryComms::BoundaryComms(const geometry::LatticeData* iLatDat,
                                 SimConfig* iSimConfig,
                                 SimulationState* iSimState,
                                 util::UnitConverter* iUnits) :
      net::IteratedAction(), mState(iSimState), mSimConfig(iSimConfig), mUnits(iUnits)
    {
      // Work out stuff for simulation (ie. should give same result on all procs)
      proc_t BCrank = 0;

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
        BCrank = topology::NetworkTopology::Instance()->GetLocalRank();

      // Since only one proc will update BCrank, the sum of all BCrank is the BCproc
      MPI_Allreduce(&BCrank, &BCproc, 1, hemelb::MpiDataType(BCrank), MPI_SUM, MPI_COMM_WORLD);

      nTotInlets = (int) mSimConfig->Inlets.size();
      nTotOutlets = (int) mSimConfig->Outlets.size();

      inletCommAlreadyChecked = new bool[nTotInlets];
      outletCommAlreadyChecked = new bool[nTotOutlets];

      for (int i = 0; i < nTotInlets; i++)
        inletCommAlreadyChecked[i] = false;
      for (int i = 0; i < nTotOutlets; i++)
        outletCommAlreadyChecked[i] = false;

      ReadParameters();

      // Work out which and how many inlets/outlets on this process

      nInlets = 0;
      nOutlets = 0;

      inlets = std::vector<int>(0);
      outlets = std::vector<int>(0);

      // Put all in/outlets onto BCproc
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        nInlets = nTotInlets;
        nOutlets = nTotOutlets;

        nInletProcs = new int[nTotInlets];
        nOutletProcs = new int[nTotOutlets];
        inletProcsList = new int*[nTotInlets];
        outletProcsList = new int*[nTotOutlets];

        for (int i = 0; i < nTotInlets; i++)
          nInletProcs[i] = 0;
        for (int i = 0; i < nTotOutlets; i++)
          nOutletProcs[i] = 0;
      }

      for (site_t i = 0; i < iLatDat->GetLocalFluidSiteCount(); i++)
      {
        if (iLatDat->GetSiteType(i) == geometry::LatticeData::INLET_TYPE
            && !util::VectorFunctions::member(inlets, iLatDat->GetBoundaryId(i)))
        {
          nInlets++;
          inlets.push_back(iLatDat->GetBoundaryId(i));
        }
        else if (iLatDat->GetSiteType(i) == geometry::LatticeData::OUTLET_TYPE
            && !util::VectorFunctions::member(outlets, iLatDat->GetBoundaryId(i)))
        {
          nOutlets++;
          outlets.push_back(iLatDat->GetBoundaryId(i));
        }
      }

      // Now BC process must find out which process belongs to what group

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        // These should be bool, but MPI only supports MPI_INT
        // For each inlet/outlet there is an array of length equal to total number of procs.
        // Each stores true/false value. True if proc of rank equal to the index contains
        // the given inlet/outlet.
        int **inletBoolList = new int*[nTotInlets];
        int **outletBoolList = new int*[nTotOutlets];

        int nProcs = topology::NetworkTopology::Instance()->GetProcessorCount();

        for (int i = 0; i < nTotInlets; i++)
        {
          inletBoolList[i] = new int[nProcs];
        }

        for (int i = 0; i < nTotOutlets; i++)
        {
          outletBoolList[i] = new int[nProcs];
        }

        MPI_Status tempStat;

        // Non-BC procs will be sending at this point
        for (int i = 0; i < nTotInlets; i++)
        {
          for (int proc = 0; proc < nProcs; proc++)
          {
            if (proc != BCproc)
            {
              MPI_Recv(&inletBoolList[i][proc], 1, MPI_INT, proc, 100, MPI_COMM_WORLD, &tempStat);
            }
            else
              inletBoolList[i][proc] = false;
          }
        }

        for (int i = 0; i < nTotOutlets; i++)
        {
          for (int proc = 0; proc < nProcs; proc++)
          {
            if (proc != BCproc)
            {
              MPI_Recv(&outletBoolList[i][proc], 1, MPI_INT, proc, 101, MPI_COMM_WORLD, &tempStat);
            }
            else
              outletBoolList[i][proc] = false;
          }
        }

        // Now we have an array for each IOlet with true (1) at indices corresponding to
        // processes that are members of that group. We have to convert this into arrays
        // of ints which store a list of processor ranks.

        for (int i = 0; i < nTotInlets; i++)
        {
          for (int j = 0; j < nProcs; j++)
          {
            if (inletBoolList[i][j])
              nInletProcs[i]++;
          }

          inletProcsList[i] = new int[nInletProcs[i]];

          int memberIndex = 0;

          for (int j = 0; j < nProcs; j++)
          {
            if (inletBoolList[i][j])
            {
              inletProcsList[i][memberIndex] = j;
              memberIndex++;
            }
          }
        }

        for (int i = 0; i < nTotOutlets; i++)
        {
          for (int j = 0; j < nProcs; j++)
          {
            if (outletBoolList[i][j])
              nOutletProcs[i]++;
          }

          outletProcsList[i] = new int[nOutletProcs[i]];

          int memberIndex = 0;

          for (int j = 0; j < nProcs; j++)
          {
            if (outletBoolList[i][j])
            {
              outletProcsList[i][memberIndex] = j;
              memberIndex++;
            }
          }
        }

        // Clear up
        for (int i = 0; i < nTotInlets; i++)
        {
          delete[] inletBoolList[i];
        }

        for (int i = 0; i < nTotOutlets; i++)
        {
          delete[] outletBoolList[i];
        }

        delete[] inletBoolList;
        delete[] outletBoolList;

      }
      else
      {
        // This is where the info about whether a proc contains a given inlet/outlet is sent
        // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
        for (int i = 0; i < nTotInlets; i++)
        {
          int inletOnThisProc = util::VectorFunctions::member(inlets, i); // true if inlet i is on this proc

          MPI_Ssend(&inletOnThisProc, 1, MPI_INT, BCproc, 100, MPI_COMM_WORLD);
        }
        for (int i = 0; i < nTotOutlets; i++)
        {
          int outletOnThisProc = util::VectorFunctions::member(outlets, i); // true if outlet i is on this proc

          MPI_Ssend(&outletOnThisProc, 1, MPI_INT, BCproc, 101, MPI_COMM_WORLD);
        }
      }

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        int nInletRequests = 0;
        int nOutletRequests = 0;

        inletRequestOffset = new int[nTotInlets];
        outletRequestOffset = new int[nTotOutlets];

        for (int i = 0; i < nTotInlets; i++)
        {
          inletRequestOffset[i] = nInletRequests;
          nInletRequests += nInletProcs[i];
        }
        for (int i = 0; i < nTotOutlets; i++)
        {
          outletRequestOffset[i] = nOutletRequests;
          nOutletRequests += nOutletProcs[i];
        }

        inlet_request = new MPI_Request[nInletRequests];
        outlet_request = new MPI_Request[nOutletRequests];
        inlet_status = new MPI_Status[nInletRequests];
        outlet_status = new MPI_Status[nOutletRequests];
      }
      else
      {
        inlet_request = new MPI_Request[nTotInlets];
        outlet_request = new MPI_Request[nTotOutlets];
        inlet_status = new MPI_Status[nTotInlets];
        outlet_status = new MPI_Status[nTotOutlets];
      }

      InitialiseBoundaryDensities();

    }

    BoundaryComms::~BoundaryComms()
    {

      // The density arrays
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        delete[] inlet_density_cycle;
        delete[] outlet_density_cycle;

        for (int i = 0; i < nTotInlets; i++)
        {
          delete[] inletProcsList[i];
          delete[] outletProcsList[i];
        }

        delete[] inletRequestOffset;
        delete[] outletRequestOffset;
        delete[] inletProcsList;
        delete[] outletProcsList;
        delete[] nInletProcs;
        delete[] nOutletProcs;
      }
      else
      {
        // On the BCproc these just point to the appropriate part of *_density_cycle array which
        // by now doesn't exist anymore
        delete[] inlet_density;
        delete[] outlet_density;
      }

      delete[] inlet_density_avg;
      delete[] outlet_density_avg;
      delete[] inlet_density_amp;
      delete[] outlet_density_amp;
      delete[] inlet_density_phs;
      delete[] outlet_density_phs;

      // Communicators and groups
      delete[] outlet_request;
      delete[] inlet_request;
      delete[] outlet_status;
      delete[] inlet_status;
    }

    void BoundaryComms::WaitForComms(const int index, unsigned int IOtype)
    {
      if (IOtype == INLET)
      {
        if (!inletCommAlreadyChecked[index])
        {
          if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
          {
            MPI_Waitall(nInletProcs[index],
                        &inlet_request[inletRequestOffset[index]],
                        &inlet_status[inletRequestOffset[index]]);
          }
          else
          {
            MPI_Wait(&inlet_request[index], &inlet_status[index]);
          }
        }
      }
      else if (IOtype == OUTLET)
      {
        if (!outletCommAlreadyChecked[index])
        {
          if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
          {
            MPI_Waitall(nOutletProcs[index],
                        &outlet_request[outletRequestOffset[index]],
                        &outlet_status[outletRequestOffset[index]]);
          }
          else
          {
            MPI_Wait(&outlet_request[index], &outlet_status[index]);
          }
        }
      }
    }

    void BoundaryComms::InitialiseBoundaryDensities()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        double w = 2.0 * PI / (double) mState->GetTimeStepsPerCycle();

        for (unsigned long time_step = 0; time_step < mState->GetTimeStepsPerCycle(); time_step++)
        {
          for (int i = 0; i < nTotInlets; i++)
          {
            inlet_density_cycle[time_step * nTotInlets + i] = inlet_density_avg[i]
                + inlet_density_amp[i] * cos(w * (double) time_step + inlet_density_phs[i]);
          }
          for (int i = 0; i < nTotOutlets; i++)
          {
            outlet_density_cycle[time_step * nTotOutlets + i] = outlet_density_avg[i]
                + outlet_density_amp[i] * cos(w * (double) time_step + outlet_density_phs[i]);
          }
        }

        inlet_density = inlet_density_cycle;
        outlet_density = outlet_density_cycle;
      }

      SendBoundaryDensities();
    }

    void BoundaryComms::RequestComms()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        unsigned long time_step = mState->GetTimeStep() % mState->GetTimeStepsPerCycle();
        inlet_density = &inlet_density_cycle[time_step * nTotInlets];
        outlet_density = &outlet_density_cycle[time_step * nTotOutlets];
      }

      SendBoundaryDensities();
    }

    void BoundaryComms::EndIteration()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        for (int i = 0; i < nTotInlets; i++)
          inletCommAlreadyChecked[i] = false;
        for (int i = 0; i < nTotOutlets; i++)
          outletCommAlreadyChecked[i] = false;
      }
      else
      {
        for (int i = 0; i < nInlets; i++)
          inletCommAlreadyChecked[inlets[i]] = false;
        for (int i = 0; i < nOutlets; i++)
          outletCommAlreadyChecked[outlets[i]] = false;
      }
    }

    void BoundaryComms::SendBoundaryDensities()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        int message = 0;

        for (int i = 0; i < nTotInlets; i++)
        {
          for (int proc = 0; proc < nInletProcs[i]; proc++)
          {
            MPI_Isend(&inlet_density[i],
                      1,
                      hemelb::MpiDataType(inlet_density[0]),
                      inletProcsList[i][proc],
                      100,
                      MPI_COMM_WORLD,
                      &inlet_request[message++]);
          }
        }

        message = 0;

        for (int i = 0; i < nTotOutlets; i++)
        {
          for (int proc = 0; proc < nOutletProcs[i]; proc++)
          {
            MPI_Isend(&outlet_density[i],
                      1,
                      hemelb::MpiDataType(outlet_density[0]),
                      outletProcsList[i][proc],
                      101,
                      MPI_COMM_WORLD,
                      &outlet_request[message++]);

          }
        }
      }
      else
      {
        for (int i = 0; i < nInlets; i++)
        {
          MPI_Irecv(&inlet_density[inlets[i]],
                    1,
                    hemelb::MpiDataType(inlet_density[0]),
                    BCproc,
                    100,
                    MPI_COMM_WORLD,
                    &inlet_request[inlets[i]]);
        }

        for (int i = 0; i < nOutlets; i++)
        {
          MPI_Irecv(&outlet_density[outlets[i]],
                    1,
                    hemelb::MpiDataType(outlet_density[0]),
                    BCproc,
                    101,
                    MPI_COMM_WORLD,
                    &outlet_request[outlets[i]]);
        }
      }
    }

    void BoundaryComms::ReadParameters()
    {
      allocateInlets();

      for (int n = 0; n < nTotInlets; n++)
      {
        hemelb::SimConfig::InOutLet *lInlet = &mSimConfig->Inlets[n];

        inlet_density_avg[n] = mUnits->ConvertPressureToLatticeUnits(lInlet->PMean) / Cs2;
        inlet_density_amp[n] = mUnits->ConvertPressureGradToLatticeUnits(lInlet->PAmp) / Cs2;
        inlet_density_phs[n] = lInlet->PPhase * DEG_TO_RAD;
      }

      allocateOutlets();

      for (int n = 0; n < nTotOutlets; n++)
      {
        hemelb::SimConfig::InOutLet *lOutlet = &mSimConfig->Outlets[n];
        outlet_density_avg[n] = mUnits->ConvertPressureToLatticeUnits(lOutlet->PMean) / Cs2;
        outlet_density_amp[n] = mUnits->ConvertPressureGradToLatticeUnits(lOutlet->PAmp) / Cs2;
        outlet_density_phs[n] = lOutlet->PPhase * DEG_TO_RAD;
      }

    }

    void BoundaryComms::allocateInlets()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        inlet_density_cycle = new distribn_t[util::NumericalFunctions::max<int>(1, nTotInlets)
            * mState->GetTimeStepsPerCycle()];
      }
      else
      {
        inlet_density = new distribn_t[nTotInlets];
      }

      inlet_density_avg = new distribn_t[nTotInlets];
      inlet_density_amp = new distribn_t[nTotInlets];
      inlet_density_phs = new distribn_t[nTotInlets];
    }

    void BoundaryComms::allocateOutlets()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        outlet_density_cycle = new distribn_t[util::NumericalFunctions::max<int>(1, nTotOutlets)
            * mState->GetTimeStepsPerCycle()];
      }
      else
      {
        outlet_density = new distribn_t[nTotOutlets];
      }

      outlet_density_avg = new distribn_t[nTotOutlets];
      outlet_density_amp = new distribn_t[nTotOutlets];
      outlet_density_phs = new distribn_t[nTotOutlets];
    }

    void BoundaryComms::Reset()
    {
      int i;

      for (i = 0; i < nTotInlets; i++)
      {
        inlet_density_avg[i] = mUnits->ConvertPressureToPhysicalUnits(inlet_density_avg[i] * Cs2);
        inlet_density_amp[i] = mUnits->ConvertPressureGradToPhysicalUnits(inlet_density_amp[i]
            * Cs2);
      }
      for (i = 0; i < nTotOutlets; i++)
      {
        outlet_density_avg[i] = mUnits->ConvertPressureToPhysicalUnits(outlet_density_avg[i] * Cs2);
        outlet_density_amp[i] = mUnits->ConvertPressureGradToPhysicalUnits(outlet_density_amp[i]
            * Cs2);
      }

      mState->DoubleTimeResolution();

      for (i = 0; i < nTotInlets; i++)
      {
        inlet_density_avg[i] = mUnits->ConvertPressureToLatticeUnits(inlet_density_avg[i]) / Cs2;
        inlet_density_amp[i] = mUnits->ConvertPressureGradToLatticeUnits(inlet_density_amp[i])
            / Cs2;
      }
      for (i = 0; i < nTotOutlets; i++)
      {
        outlet_density_avg[i] = mUnits->ConvertPressureToLatticeUnits(outlet_density_avg[i]) / Cs2;
        outlet_density_amp[i] = mUnits->ConvertPressureGradToLatticeUnits(outlet_density_amp[i])
            / Cs2;
      }

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        delete[] inlet_density_cycle;
        delete[] outlet_density_cycle;

        inlet_density_cycle = new distribn_t[hemelb::util::NumericalFunctions::max<int>(1,
                                                                                        nTotInlets)
            * mState->GetTimeStepsPerCycle()];
        outlet_density_cycle
            = new distribn_t[hemelb::util::NumericalFunctions::max<int>(1, nTotOutlets)
                * mState->GetTimeStepsPerCycle()];
      }

      InitialiseBoundaryDensities();
    }

    // Calculate the BCs for each boundary site type and the
    // non-equilibrium distribution functions.
    void BoundaryComms::CalculateBC(distribn_t f[],
                                    hemelb::geometry::LatticeData::SiteType iSiteType,
                                    unsigned int iBoundaryId,
                                    distribn_t *density,
                                    distribn_t *vx,
                                    distribn_t *vy,
                                    distribn_t *vz,
                                    distribn_t f_neq[])
    {
      distribn_t dummy_density;

      for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
      {
        f_neq[l] = f[l];
      }

      if (iSiteType == hemelb::geometry::LatticeData::FLUID_TYPE)
      {
        D3Q15::CalculateDensityAndVelocity(f, *density, *vx, *vy, *vz);
      }
      else
      {
        if (iSiteType == hemelb::geometry::LatticeData::INLET_TYPE)
        {
          *density = inlet_density[iBoundaryId];
        }
        else
        {
          *density = outlet_density[iBoundaryId];
        }

        D3Q15::CalculateDensityAndVelocity(f, dummy_density, *vx, *vy, *vz);
        D3Q15::CalculateFeq(*density, *vx, *vy, *vz, f);

      }
      for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
      {
        f_neq[l] -= f[l];
      }

    }

  }
}
