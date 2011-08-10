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

      ReadParameters();

      inlet_groups = new MPI_Group[nTotInlets];
      inlet_comms = new MPI_Comm[nTotInlets];
      outlet_groups = new MPI_Group[nTotOutlets];
      outlet_comms = new MPI_Comm[nTotOutlets];

      // Work out which and how many inlets/outlets on this process

      nInlets = 0;
      nOutlets = 0;

      inlets = std::vector<int>(0);
      outlets = std::vector<int>(0);

      // Put all in/outlets onto BCproc
      if (topology::NetworkTopology::Instance()->GetLocalRank() == BCproc)
      {
        nInlets = nTotInlets;
        for (int i = 0; i < nTotInlets; i++)
          inlets.push_back((int) i);

        nOutlets = nTotOutlets;
        for (int i = 0; i < nTotOutlets; i++)
          outlets.push_back((int) i);
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

      // Have to be sorted in order to prevent deadlocks when broadcasting
      // e.g. One process has two inlets, say 3 and 7. BCproc will always have them in order
      // so it will attempt to broadcast inlet 3 data first. If the process with 3 and 7 asks
      // for inlet 7 first we will have a deadlock, since it will be waiting for a broadcast
      // on a different communicator.
      util::VectorFunctions::BubbleSort(inlets);
      util::VectorFunctions::BubbleSort(outlets);

      // Now all processes must find out which process belongs to what group

      // These should be bool, but MPI only supports MPI_INT
      // For each inlet/outlet there is an array of length equal to total number of procs.
      // Each stores true/false value. True if proc of rank equal to the index contains
      // the given inlet/outlet.
      int **inletProcsList = new int*[nTotInlets];
      int **outletProcsList = new int*[nTotOutlets];

      int nProcs = topology::NetworkTopology::Instance()->GetProcessorCount();

      for (int i = 0; i < nTotInlets; i++)
      {
        inletProcsList[i] = new int[nProcs];
      }

      for (int i = 0; i < nTotOutlets; i++)
      {
        outletProcsList[i] = new int[nProcs];
      }

      // This is where the info about whether a proc contains a given inlet/outlet is sent/received
      // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
      for (int i = 0; i < nTotInlets; i++)
      {
        int inletOnThisProc = util::VectorFunctions::member(inlets, i); // true if inlet i is on this proc

        MPI_Allgather(&inletOnThisProc, 1, MPI_INT, inletProcsList[i], 1, MPI_INT, MPI_COMM_WORLD);
      }
      for (int i = 0; i < nTotOutlets; i++)
      {
        int outletOnThisProc = util::VectorFunctions::member(outlets, i); // true if outlet i is on this proc

        MPI_Allgather(&outletOnThisProc, 1, MPI_INT, outletProcsList[i], 1, MPI_INT, MPI_COMM_WORLD);
      }

      // Now we have an array for each group with true (1) at indices corresponding to
      // processes that are members of that group. We have to convert this into arrays
      // of ints which store a list of processor ranks.

      int **inletGroupMembers = new int*[nTotInlets];
      int **outletGroupMembers = new int*[nTotOutlets];

      MPI_Group orig_group;

      MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

      for (int i = 0; i < nTotInlets; i++)
      {
        int nGroupMembers = 0;

        for (int j = 0; j < nProcs; j++)
        {
          if (inletProcsList[i][j])
            nGroupMembers++;
        }

        inletGroupMembers[i] = new int[nGroupMembers];

        int memberIndex = 0;

        for (int j = 0; j < nProcs; j++)
        {
          if (inletProcsList[i][j])
          {
            inletGroupMembers[i][memberIndex] = j;
            memberIndex++;
          }
        }

        // Create the group and comm now
        MPI_Group_incl(orig_group, nGroupMembers, inletGroupMembers[i], &inlet_groups[i]);
        MPI_Comm_create(MPI_COMM_WORLD, inlet_groups[i], &inlet_comms[i]);
      }

      for (int i = 0; i < nTotOutlets; i++)
      {
        int nGroupMembers = 0;

        for (int j = 0; j < nProcs; j++)
        {
          if (outletProcsList[i][j])
            nGroupMembers++;
        }

        outletGroupMembers[i] = new int[nGroupMembers];

        int memberIndex = 0;

        for (int j = 0; j < nProcs; j++)
        {
          if (outletProcsList[i][j])
          {
            outletGroupMembers[i][memberIndex] = j;
            memberIndex++;
          }
        }

        // Create the group and comm now
        MPI_Group_incl(orig_group, nGroupMembers, outletGroupMembers[i], &outlet_groups[i]);
        MPI_Comm_create(MPI_COMM_WORLD, outlet_groups[i], &outlet_comms[i]);
      }

      InitialiseBoundaryDensities();

      // Clear up

      for (int i = 0; i < nTotInlets; i++)
      {
        delete[] inletProcsList[i];
        delete[] inletGroupMembers[i];
      }

      for (int i = 0; i < nTotOutlets; i++)
      {
        delete[] outletProcsList[i];
        delete[] outletGroupMembers[i];
      }

      delete[] inletProcsList;
      delete[] inletGroupMembers;
      delete[] outletProcsList;
      delete[] outletGroupMembers;

    }

    BoundaryComms::~BoundaryComms()
    {

      // The density arrays
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        delete[] inlet_density_cycle;
        delete[] outlet_density_cycle;
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
      delete[] outlet_comms;
      delete[] inlet_comms;
      delete[] outlet_groups;
      delete[] inlet_groups;
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

      BroadcastBoundaryDensities();
    }

    void BoundaryComms::RequestComms()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        unsigned long time_step = mState->GetTimeStep() % mState->GetTimeStepsPerCycle();
        inlet_density = &inlet_density_cycle[time_step * nTotInlets];
        outlet_density = &outlet_density_cycle[time_step * nTotOutlets];
      }
      BroadcastBoundaryDensities();
    }

    void BoundaryComms::BroadcastBoundaryDensities()
    {
      for (int i = 0; i < nInlets; i++)
      {
        MPI_Bcast(&inlet_density[inlets[i]],
                  1,
                  hemelb::MpiDataType(inlet_density[0]),
                  BCproc,
                  inlet_comms[inlets[i]]);
      }
      for (int i = 0; i < nOutlets; i++)
      {
        MPI_Bcast(&outlet_density[outlets[i]],
                  1,
                  hemelb::MpiDataType(outlet_density[0]),
                  BCproc,
                  outlet_comms[outlets[i]]);
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
