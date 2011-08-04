#ifndef HEMELB_LB_BOUNDARYCOMMS_H
#define HEMELB_LB_BOUNDARYCOMMS_H

#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {

    class BoundaryComms
    {
      public:
        BoundaryComms(const geometry::LatticeData* iLatDat, const SimConfig* iSimConfig);
        ~BoundaryComms();

        void UpdateBoundaryDensities(unsigned long time_step,
                                     unsigned long timeStepsPerCycle,
                                     distribn_t* inlet_density,
                                     distribn_t* inlet_density_avg,
                                     distribn_t* inlet_density_amp,
                                     distribn_t* inlet_density_phs,
                                     distribn_t* outlet_density,
                                     distribn_t* outlet_density_avg,
                                     distribn_t* outlet_density_amp,
                                     distribn_t* outlet_density_phs);

        void BroadcastBoundaryDensities(distribn_t* inlet_density, distribn_t* outlet_density);

      private:
        proc_t BCproc; // Process responsible for sending out BC info

        // Total number of inlets/outlets in simulation
        size_t nTotInlets;
        size_t nTotOutlets;

        // Number of inlets/outlets on this process
        size_t nInlets;
        size_t nOutlets;

        // List of indices of inlets/outlets on this process
        std::vector<int>* inlets;
        std::vector<int>* outlets;

        // Communicators and groups
        MPI_Comm* outlet_comms;
        MPI_Comm* inlet_comms;
        MPI_Group* outlet_groups;
        MPI_Group* inlet_groups;

        template<typename T>
        bool member(std::vector<T> &list, T element);

        // Just for testing ATM
        void printStuff();

    };

    template<typename T>
    bool BoundaryComms::member(std::vector<T> &list, T element)
    {
      for (unsigned int i = 0; i < list.size(); i++)
      {
        if (element == list[i])
          return true;
      }

      return false;
    }

  }
}

#endif /* HEMELB_LB_BOUNDARYCOMMS_H */
