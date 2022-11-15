// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LB_H
#define HEMELB_LB_LB_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "net/IOCommunicator.h"
#include "lb/SimulationState.h"
#include "lb/InitialCondition.h"
#include "lb/iolets/BoundaryValues.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/UnitConverter.h"
#include "reporting/Timers.h"
#include "lb/BuildSystemInterface.h"
#include "Traits.h"
#include <typeinfo>

namespace hemelb
{
  /**
   * Namespace 'lb' contains classes for the scientific core of the Lattice Boltzman simulation
   */
  namespace lb
  {
    /**
     * Class providing core Lattice Boltzmann functionality.
     * Implements the IteratedAction interface.
     */
    template<class TRAITS = hemelb::Traits<>>
    class LBM : public net::IteratedAction
    {
      public:
        //! Instanciation type
        typedef TRAITS Traits;
      private:
        typedef typename Traits::Lattice LatticeType;
        // Use the kernel specified through the build system. This will select one of the above classes.
        typedef typename Traits::Kernel LBKernel;

        typedef typename Traits::Streamer tMidFluidCollision;
        // Use the wall boundary condition specified through the build system.
        typedef typename Traits::WallBoundary tWallCollision;
        // Use the inlet BC specified by the build system
        typedef typename Traits::InletBoundary tInletCollision;
        // Use the outlet BC specified by the build system
        typedef typename Traits::OutletBoundary tOutletCollision;
        // And again but for sites that are both in-/outlet and wall
        typedef typename Traits::WallInletBoundary tInletWallCollision;
        typedef typename Traits::WallOutletBoundary tOutletWallCollision;

      public:
        /**
         * Constructor, stage 1.
         * Object so initialized is not ready for simulation.
         * Must have Initialise(...) called also. Constructor separated due to need to access
         * the partially initialized LBM in order to initialize the arguments to the second construction phase.
         */
        LBM(LbmParameters params, net::Net* net,
            geometry::FieldData* latDat, SimulationState* simState, reporting::Timers &atimings,
            geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager);
        ~LBM();

        void RequestComms(); ///< part of IteratedAction interface.
        void PreSend(); ///< part of IteratedAction interface.
        void PreReceive(); ///< part of IteratedAction interface.
        void PostReceive(); ///< part of IteratedAction interface.
        void EndIteration(); ///< part of IteratedAction interface.

        site_t TotalFluidSiteCount() const;
        void SetTotalFluidSiteCount(site_t);

        /**
         * Second constructor.
         *
         */
        void Initialise(iolets::BoundaryValues* iInletValues,
                        iolets::BoundaryValues* iOutletValues,
                        const util::UnitConverter* iUnits);

        void SetInitialConditions(lb::InitialCondition const& ic_conf, const net::IOCommunicator& ioComms);

        hemelb::lb::LbmParameters *GetLbmParams();
        lb::MacroscopicPropertyCache& GetPropertyCache();

      private:

        void InitCollisions();
        // The following function pair simplify initialising the site ranges for each collider object.
        void InitInitParamsSiteRanges(kernels::InitParams& initParams, unsigned& state);
        void AdvanceInitParamsSiteRanges(kernels::InitParams& initParams, unsigned& state);
        /**
         * Ensure that the BoundaryValues objects have all necessary fields populated.
         */
        void PrepareBoundaryObjects();

        void ReadParameters();

        void handleIOError(int iError);

        // Collision objects
        tMidFluidCollision* mMidFluidCollision;
        tWallCollision* mWallCollision;
        tInletCollision* mInletCollision;
        tOutletCollision* mOutletCollision;
        tInletWallCollision* mInletWallCollision;
        tOutletWallCollision* mOutletWallCollision;

        template<typename Collision>
        void StreamAndCollide(Collision* collision, const site_t iFirstIndex,
                              const site_t iSiteCount)
        {
            collision->StreamAndCollide(iFirstIndex,
                                                 iSiteCount,
                                                 &mParams,
                                                 *mLatDat,
                                                 propertyCache);
        }

        template<typename Collision>
        void PostStep(Collision* collision, const site_t iFirstIndex, const site_t iSiteCount)
        {
            collision->DoPostStep(iFirstIndex,
                                           iSiteCount,
                                           &mParams,
                                           *mLatDat,
                                           propertyCache);

        }

        net::Net* mNet;
        geometry::FieldData* mLatDat;
        SimulationState* mState;
        iolets::BoundaryValues *mInletValues, *mOutletValues;

        LbmParameters mParams;

        const util::UnitConverter* mUnits;

        hemelb::reporting::Timers &timings;

        MacroscopicPropertyCache propertyCache;

        geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager;
    };

  } // Namespace lb
} // Namespace hemelb
#endif // HEMELB_LB_LB_H
