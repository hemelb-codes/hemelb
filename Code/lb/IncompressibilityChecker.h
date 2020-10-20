
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_INCOMPRESSIBILITYCHECKER_H
#define HEMELB_LB_INCOMPRESSIBILITYCHECKER_H
#include "comm/Communicator.h"
#include "timestep/CollectiveActor.h"
#include "geometry/LatticeData.h"
#include "lb/MacroscopicPropertyCache.h"
#include "reporting/Reportable.h"

namespace hemelb
{
  namespace lb
  {
    /**
     * Domain reference density used to compute relative differences.
     *
     * @todo #23 This value should not be hardcoded here but read from UnitConverter
     */
    static const distribn_t REFERENCE_DENSITY = 1.0;

    struct DensityEtc {
        double min;
        double max;
        double maxVel;
    };

    class IncompressibilityChecker : public timestep::CollectiveActor,
                                     public reporting::Reportable
    {
      public:
        /**
         * Constructor
         *
         * @param latticeData geometry object
         * @param comms communicator object
         * @param maximumRelativeDensityDifferenceAllowed maximum density difference allowed in the domain (relative to reference density, default 5%)
         */
      IncompressibilityChecker(const geometry::LatticeData * latticeData, comm::Communicator::ConstPtr comms,
                                 lb::MacroscopicPropertyCache& propertyCache,
                                 reporting::Timers& timings,
                                 distribn_t maximumRelativeDensityDifferenceAllowed = 0.05);

        /**
         * Destructor
         */
        virtual ~IncompressibilityChecker();

        void Report(reporting::Dict& dictionary) override;

        /**
         * Returns smallest density in the domain as agreed by all the processes.
         *
         * @return smallest density
         */
        distribn_t GetGlobalSmallestDensity() const;

        /**
         * Returns largest density in the domain as agreed by all the processes.
         *
         * @return largest density
         */
        distribn_t GetGlobalLargestDensity() const;

        /**
         * Return current maximum density difference across the domain (relative to domain reference density).
         *
         * @return current maximum density difference across the domain
         */
        double GetMaxRelativeDensityDifference() const;

        /**
         * Return allowed maximum density difference across the domain (relative to domain reference density).
         *
         * @return allowed maximum density difference across the domain
         */
        double GetMaxRelativeDensityDifferenceAllowed() const;

        /**
         * Checks whether the maximum density difference is smaller that the maximum allowed.
         *
         * @return whether the maximum density difference is smaller that the maximum allowed
         */
        bool IsDensityDiffWithinRange() const;

        /**
         * Return largest velocity magnitude in the domain as agreed by all the processes.
         *
         * @return largest velocity magnitude
         */
        double GetGlobalLargestVelocityMagnitude() const;

        static MPI_Datatype GetDensityMpiDatatype();

      private:

      // (Collective)Actor interface
      inline void BeginAll() override {}
      inline void Begin() override {}
      // Receive provided by CollectiveActor
      // Compute the local state
      void PreSend() override;
      // Initiate the collective
      void Send() override;
      inline void PreWait() override {}
      // Wait provided by CollectiveActor
      inline void End() override {}
      inline void EndAll() override {}

        static void UpdateDensityEtc(const DensityEtc& in, DensityEtc& inout);
        static void MpiOpUpdateFunc(void* invec, void* inoutvec,
                                    int *len, MPI_Datatype *datatype);
        /** Pointer to lattice data object. */
        const geometry::LatticeData * mLatDat;

        /** Cache of macroscopic properties (including density). */
        lb::MacroscopicPropertyCache& propertyCache;

        /** Pointer to the simulation state used in the rest of the simulation. */
        // lb::SimulationState* mSimState;

        /** Maximum density difference allowed in the domain (relative to reference density) */
        distribn_t maximumRelativeDensityDifferenceAllowed;

        /** Density tracker with the densities agreed on. */
        DensityEtc globalDensity;

        /** Density tracker of this node and its children to propagate upwards. */
        DensityEtc localDensity;

        /** Custom operator for reduction/ */
        MPI_Op reduction;

        /** Time spent checking stuff */
        reporting::Timer& workTimer;
    };

  }
}

#endif /* HEMELB_LB_INCOMPRESSIBILITYCHECKER_H */
