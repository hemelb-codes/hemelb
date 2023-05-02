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

/**
 * Namespace 'lb' contains classes for the scientific core of the Lattice Boltzmann simulation
 */
namespace hemelb::lb
{
    /**
     * Class providing core Lattice Boltzmann functionality.
     * Implements the IteratedAction interface.
     */
    template<class TRAITS = hemelb::Traits<>>
    class LBM : public net::IteratedAction
    {
      public:
        //! Instantiation type
        using Traits = TRAITS;
      private:
        using LatticeType = typename Traits::Lattice;
        // Use the kernel specified through the build system. This will select one of the above classes.
        using LBKernel = typename Traits::Kernel;

        using tMidFluidCollision = typename Traits::Streamer;
        // Use the wall boundary condition specified through the build system.
        using tWallCollision = typename Traits::WallBoundary;
        // Use the inlet BC specified by the build system
        using tInletCollision = typename Traits::InletBoundary;
        // Use the outlet BC specified by the build system
        using tOutletCollision = typename Traits::OutletBoundary;
        // And again but for sites that are both in-/outlet and wall
        using tInletWallCollision = typename Traits::WallInletBoundary;
        using tOutletWallCollision = typename Traits::WallOutletBoundary;

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
        ~LBM() override = default;

        void RequestComms() override; ///< part of IteratedAction interface.
        void PreSend() override; ///< part of IteratedAction interface.
        void PreReceive() override; ///< part of IteratedAction interface.
        void PostReceive() override; ///< part of IteratedAction interface.
        void EndIteration() override; ///< part of IteratedAction interface.

        [[nodiscard]] site_t TotalFluidSiteCount() const;
        void SetTotalFluidSiteCount(site_t);

        /**
         * Second constructor.
         *
         */
        void Initialise(iolets::BoundaryValues* iInletValues,
                        iolets::BoundaryValues* iOutletValues);

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
        // TODO: these should probably be optional but we can't because some of the collisions aren't copyable
        std::unique_ptr<tMidFluidCollision> mMidFluidCollision;
        std::unique_ptr<tWallCollision> mWallCollision;
        std::unique_ptr<tInletCollision> mInletCollision;
        std::unique_ptr<tOutletCollision> mOutletCollision;
        std::unique_ptr<tInletWallCollision> mInletWallCollision;
        std::unique_ptr<tOutletWallCollision> mOutletWallCollision;

        template<typename Collision>
        void StreamAndCollide(Collision& collision, const site_t iFirstIndex,
                              const site_t iSiteCount)
        {
            collision.StreamAndCollide(iFirstIndex,
                                                 iSiteCount,
                                                 &mParams,
                                                 *mLatDat,
                                                 propertyCache);
        }

        template<typename Collision>
        void PostStep(Collision& collision, const site_t iFirstIndex, const site_t iSiteCount)
        {
            collision.DoPostStep(iFirstIndex,
                                           iSiteCount,
                                           &mParams,
                                           *mLatDat,
                                           propertyCache);

        }

        net::Net* mNet;
        geometry::FieldData* mLatDat;
        SimulationState* mState;
        iolets::BoundaryValues *mInletValues = nullptr,
                *mOutletValues = nullptr;

        LbmParameters mParams;

        hemelb::reporting::Timers &timings;

        MacroscopicPropertyCache propertyCache;

        geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager;
    };

}
#endif // HEMELB_LB_LB_H
