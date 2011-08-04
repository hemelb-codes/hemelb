#include "lb/BoundaryComms.h"
#include "topology/NetworkTopology.h"
#include <math.h>

namespace hemelb
{
  namespace lb
  {

    BoundaryComms::BoundaryComms(const geometry::LatticeData* iLatDat, const SimConfig* mSimConfig)
    {
      // Work out stuff for simulation (ie. should give same result on all procs)
      BCproc = 0;

      nTotInlets = mSimConfig->Inlets.size();
      nTotOutlets = mSimConfig->Outlets.size();

      inlet_groups = new MPI_Group[nTotInlets];
      inlet_comms = new MPI_Comm[nTotInlets];
      outlet_groups = new MPI_Group[nTotOutlets];
      outlet_comms = new MPI_Comm[nTotOutlets];

      // Work out which and how many inlets/outlets on this process

      nInlets = 0;
      nOutlets = 0;

      inlets = new std::vector<int>(0);
      outlets = new std::vector<int>(0);

      // Put all in/outlets onto BCproc
      if (topology::NetworkTopology::Instance()->GetLocalRank() == 0)
      {
        nInlets = nTotInlets;
        for (size_t i = 0; i < nTotInlets; i++)
          inlets->push_back((int) i);

        nOutlets = nTotOutlets;
        for (size_t i = 0; i < nTotOutlets; i++)
          outlets->push_back((int) i);
      }

      for (site_t i = 0; i < iLatDat->GetLocalFluidSiteCount(); i++)
      {
        if (iLatDat->GetSiteType(i) == geometry::LatticeData::INLET_TYPE
            && !member(*inlets, iLatDat->GetBoundaryId(i)))
        {
          nInlets++;
          inlets->push_back(iLatDat->GetBoundaryId(i));
        }
        else if (iLatDat->GetSiteType(i) == geometry::LatticeData::OUTLET_TYPE
            && !member(*outlets, iLatDat->GetBoundaryId(i)))
        {
          nOutlets++;
          outlets->push_back(iLatDat->GetBoundaryId(i));
        }
      }

      // Now all processes must find out which process belongs to what group

      // These should be bool, but MPI only supports MPI_INT
      // For each inlet/outlet there is an array of length equal to total number of procs.
      // Each stores true/false value. True if proc of rank equal to the index contains
      // the given inlet/outlet.
      int **inletProcsList = new int*[nTotInlets];
      int **outletProcsList = new int*[nTotOutlets];

      int nProcs = topology::NetworkTopology::Instance()->GetProcessorCount();

      for (size_t i = 0; i < nTotInlets; i++)
      {
        inletProcsList[i] = new int[nProcs];
      }

      for (size_t i = 0; i < nTotOutlets; i++)
      {
        outletProcsList[i] = new int[nProcs];
      }

      // This is where the info about whether a proc contains a given inlet/outlet is sent/received
      // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
      for (size_t i = 0; i < nTotInlets; i++)
      {
        int inletOnThisProc = member<int> (*inlets, (int) i); // true if inlet i is on this proc

        MPI_Allgather(&inletOnThisProc, 1, MPI_INT, inletProcsList[i], 1, MPI_INT, MPI_COMM_WORLD);
      }
      for (size_t i = 0; i < nTotOutlets; i++)
      {
        int outletOnThisProc = member<int> (*outlets, (int) i); // true if inlet i is on this proc

        MPI_Allgather(&outletOnThisProc, 1, MPI_INT, outletProcsList[i], 1, MPI_INT, MPI_COMM_WORLD);
      }

      // Now we have an array for each group with true (1) at indices corresponding to
      // processes that are members of that group. We have to convert this into arrays
      // of ints which store a list of processor ranks.

      int **inletGroupMembers = new int*[nTotInlets];
      int **outletGroupMembers = new int*[nTotOutlets];

      MPI_Group orig_group;

      MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

      for (size_t i = 0; i < nTotInlets; i++)
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

        // Create the group now
        MPI_Group_incl(orig_group, nGroupMembers, inletGroupMembers[i], &inlet_groups[i]);
      }

      for (size_t i = 0; i < nTotOutlets; i++)
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

        // Create the group now
        MPI_Group_incl(orig_group, nGroupMembers, outletGroupMembers[i], &outlet_groups[i]);
      }

      // Finally create the comms
      for (size_t i = 0; i < nTotInlets; i++)
      {
        MPI_Comm_create(MPI_COMM_WORLD, inlet_groups[i], &inlet_comms[i]);
      }
      for (size_t i = 0; i < nTotOutlets; i++)
      {
        MPI_Comm_create(MPI_COMM_WORLD, outlet_groups[i], &outlet_comms[i]);
      }

      // Clear up

      for (size_t i = 0; i < nTotInlets; i++)
      {
        delete[] inletProcsList[i];
        delete[] inletGroupMembers[i];
      }

      for (size_t i = 0; i < nTotOutlets; i++)
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
      delete inlets;
      delete outlets;

      // Communicators and groups
      delete[] outlet_comms;
      delete[] inlet_comms;
      delete[] outlet_groups;
      delete[] inlet_groups;
    }

    void BoundaryComms::UpdateBoundaryDensities(unsigned long time_step,
                                                unsigned long timeStepsPerCycle,
                                                distribn_t* inlet_density,
                                                distribn_t* inlet_density_avg,
                                                distribn_t* inlet_density_amp,
                                                distribn_t* inlet_density_phs,
                                                distribn_t* outlet_density,
                                                distribn_t* outlet_density_avg,
                                                distribn_t* outlet_density_amp,
                                                distribn_t* outlet_density_phs)
    {
      double w = 2.0 * PI / (double) timeStepsPerCycle;

      for (size_t i = 0; i < nTotInlets; i++)
      {
        inlet_density[i] = inlet_density_avg[i] + inlet_density_amp[i] * cos(w * (double) time_step
            + inlet_density_phs[i]);
      }
      for (size_t i = 0; i < nTotOutlets; i++)
      {
        outlet_density[i] = outlet_density_avg[i] + outlet_density_amp[i] * cos(w
            * (double) time_step + outlet_density_phs[i]);
      }
    }

    void BoundaryComms::BroadcastBoundaryDensities(distribn_t* inlet_density,
                                                   distribn_t* outlet_density)
    {
      for (size_t i = 0; i < nInlets; i++)
      {
        MPI_Bcast(&inlet_density[ (*inlets)[i]],
                  1,
                  hemelb::MpiDataType(inlet_density[0]),
                  BCproc,
                  inlet_comms[ (*inlets)[i]]);
      }
      for (size_t i = 0; i < nOutlets; i++)
      {
        MPI_Bcast(&outlet_density[ (*outlets)[i]],
                  1,
                  hemelb::MpiDataType(outlet_density[0]),
                  BCproc,
                  outlet_comms[ (*outlets)[i]]);
      }
    }

    void BoundaryComms::printStuff()
    {
      std::cout << topology::NetworkTopology::Instance()->GetLocalRank() << " nInlets: " << nInlets
          << "/" << nTotInlets << " nOutlets: " << nOutlets << "/" << nTotOutlets << " inlets: ";

      for (unsigned int i = 0; i < inlets->size(); i++)
      {
        std::cout << (*inlets)[i] << " ";
      }

      std::cout << "outlets: ";

      for (unsigned int i = 0; i < outlets->size(); i++)
      {
        std::cout << (*outlets)[i] << " ";
      }

      std::cout << std::endl;
    }

  }
}
