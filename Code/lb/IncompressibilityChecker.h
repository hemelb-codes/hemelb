
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_INCOMPRESSIBILITYCHECKER_H
#define HEMELB_LB_INCOMPRESSIBILITYCHECKER_H

#include "geometry/LatticeData.h"
#include "lb/MacroscopicPropertyCache.h"
#include "net/PhasedBroadcastRegular.h"
#include "reporting/Reportable.h"
#include <cfloat>

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

    template<class BroadcastPolicy>
    class IncompressibilityChecker : public BroadcastPolicy,
                                     public reporting::Reportable
    {
        /**
         * This class uses the phased broadcast infrastructure to keep track of the maximum density difference across the domain.
         */

      public:
        class DensityTracker
        {
            /**
             * This is convenience class that encapsulates fluid magnitudes of interest being tracked across the domain. At the
             * moment it tracks maximum and minimum density, it can be extended to accommodate other magnitudes.
             */
          public:
            /** Number of densities being tracked. Needs to be exposed for send/receive. */
            static const unsigned DENSITY_TRACKER_SIZE = 3;

            /** Identifiers of the densities being tracked. Cardinality must be kept consistent with DENSITY_TRACKER_SIZE */
            typedef enum
            {
              MIN_DENSITY = 0u,
              MAX_DENSITY,
              MAX_VELOCITY_MAGNITUDE
            } DensityTrackerIndices;

            /**
             * Default constructor, initialises the tracker with some large/small min/max densities.
             */
            DensityTracker();

            /**
             * Constructor, initialises the tracker with density values provided.
             *
             * @param densityValues density values used for initialisation
             */
            DensityTracker(distribn_t* const densityValues);

            /**
             * Destructor
             */
            ~DensityTracker();

            /**
             * Updates densities based on the values of another DensityTracker object.
             *
             * @param newValues density tracker object to copy from
             */
            void operator=(const DensityTracker& newValues);

            /**
             * Access individually each of the densities tracked.
             *
             * @param index index of the density of interest (see DensityTrackerIndices enum for a list)
             * @return density value
             */
            distribn_t& operator[](DensityTrackerIndices densityIndex) const;

            /**
             * Returns a pointer to the internal array used to store densities.
             * Only to be used for sending/receiving. Use [] operator, instead.
             *
             * @return pointer to the internal array used to store densities
             */
            distribn_t* GetDensitiesArray() const;

            /**
             * Updates min/max densities with values in newValue object if they are smaller/larger
             *
             * @param newValues new values to be considered for an update
             */
            void UpdateDensityTracker(const DensityTracker& newValues);

            /**
             * Updates min/max densities with value newValue if it is smaller/larger
             *
             * @param newDensity new density value to be considered for an update
             * @param newVelocityMagnitude new velocity magnitude to be considered for an update
             */
            void UpdateDensityTracker(distribn_t newDensity, distribn_t newVelocityMagnitude);

          private:
            /** Array storing all the densities being tracked */
            distribn_t* densitiesArray;

            /** Keeps track on whether densitiesArray was allocated locally */
            bool allocatedHere;
        };

        /**
         * Constructor
         *
         * @param latticeData geometry object
         * @param net network interface object
         * @param simState simulation state
         * @param maximumRelativeDensityDifferenceAllowed maximum density difference allowed in the domain (relative to reference density, default 5%)
         */
        IncompressibilityChecker(const geometry::LatticeData * latticeData,
                                 net::Net* net,
                                 SimulationState* simState,
                                 lb::MacroscopicPropertyCache& propertyCache,
                                 reporting::Timers& timings,
                                 distribn_t maximumRelativeDensityDifferenceAllowed = 0.05);

        /**
         * Destructor
         */
        virtual ~IncompressibilityChecker();

        void Report(ctemplate::TemplateDictionary& dictionary);

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
         * Returns whether the first max/min density reduction operation has finished and
         * therefore there are density values available.
         *
         * @return whether there are density values available
         */
        bool AreDensitiesAvailable() const;

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

      protected:
        /**
         * Override the methods from the base class to propagate data from the root, and
         * to send data about this node and its childrens' density trackers up towards the root.
         */
        void ProgressFromChildren(unsigned long splayNumber);
        void ProgressFromParent(unsigned long splayNumber);
        void ProgressToChildren(unsigned long splayNumber);
        void ProgressToParent(unsigned long splayNumber);

        /**
         * Take the combined density tracker information and start passing it back down the tree.
         */
        void TopNodeAction();

        /**
         * Override the method from the base class to use the data from child nodes.
         */
        void PostReceiveFromChildren(unsigned long splayNumber);

        /**
         * Override the method from the base class for a node in the tree to update its density tracker
         *
         * @param splayNumber The parameter splayNumber is 0 indexed and less than splay (a template parameter of PhasedBroadcast)
         */
        virtual void PostSendToParent(unsigned long splayNumber);

        /**
         * Use the density tracker sent by the root node.
         */
        void Effect();

      private:

        /**
         * Slightly arbitrary spread factor for the tree.
         *
         * @todo #23 This is defined in StabilityChecker as well, refactor somewhere else?
         */
        static const unsigned int SPREADFACTOR = 10u;

        /** Pointer to lattice data object. */
        const geometry::LatticeData * mLatDat;

        /** Cache of macroscopic properties (including density). */
        lb::MacroscopicPropertyCache& propertyCache;

        /** Pointer to the simulation state used in the rest of the simulation. */
        lb::SimulationState* mSimState;

        /** Timing object. */
        reporting::Timers& timings;

        /** Maximum density difference allowed in the domain (relative to reference density) */
        distribn_t maximumRelativeDensityDifferenceAllowed;

        /** Density tracker with the densities agreed on. */
        DensityTracker* globalDensityTracker;

        /** Density tracker of this node and its children to propagate upwards. */
        DensityTracker upwardsDensityTracker;

        /** Density tracker as agreed at the root node, to pass downwards. */
        DensityTracker downwardsDensityTracker;

        /** Array for storing the passed-up densities from child nodes. */
        distribn_t childrenDensitiesSerialised[SPREADFACTOR * DensityTracker::DENSITY_TRACKER_SIZE];
    };

  }
}

#endif /* HEMELB_LB_INCOMPRESSIBILITYCHECKER_H */
