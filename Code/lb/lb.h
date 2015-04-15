// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_LB_H
#define HEMELB_LB_LB_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "net/IOCommunicator.h"
#include "lb/SimulationState.h"
#include "lb/iolets/BoundaryValues.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/UnitConverter.h"
#include "configuration/SimConfig.h"
#include "reporting/Timers.h"
#include "lb/BuildSystemInterface.h"
#include "traits.h"
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
        LBM(hemelb::configuration::SimConfig *iSimulationConfig, net::Net* net,
            geometry::LatticeData* latDat, SimulationState* simState, reporting::Timers &atimings,
            geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager);
        ~LBM();

        void RequestComms(); ///< part of IteratedAction interface.
        void PreSend(); ///< part of IteratedAction interface.
        void PreReceive(); ///< part of IteratedAction interface.
        void PostReceive(); ///< part of IteratedAction interface.
        void EndIteration(); ///< part of IteratedAction interface.

        site_t TotalFluidSiteCount() const;
        void SetTotalFluidSiteCount(site_t);
        int InletCount() const
        {
          return inletCount;
        }
        int OutletCount() const
        {
          return outletCount;
        }

        /**
         * Second constructor.
         *
         */
        void Initialise(vis::Control* iControl, iolets::BoundaryValues* iInletValues,
                        iolets::BoundaryValues* iOutletValues, const util::UnitConverter* iUnits);

        void ReadVisParameters();

        void CalculateMouseFlowField(const ScreenDensity densityIn, const ScreenStress stressIn,
                                     const LatticeDensity density_threshold_min,
                                     const LatticeDensity density_threshold_minmax_inv,
                                     const LatticeStress stress_threshold_max_inv,
                                     PhysicalPressure &mouse_pressure,
                                     PhysicalStress &mouse_stress);

        hemelb::lb::LbmParameters *GetLbmParams();
        lb::MacroscopicPropertyCache& GetPropertyCache();

      private:
        void SetInitialConditions();

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
          if (mVisControl->IsRendering())
          {
            collision->template StreamAndCollide<true>(iFirstIndex,
                                                       iSiteCount,
                                                       &mParams,
                                                       mLatDat,
                                                       propertyCache);
          }
          else
          {
            collision->template StreamAndCollide<false>(iFirstIndex,
                                                        iSiteCount,
                                                        &mParams,
                                                        mLatDat,
                                                        propertyCache);
          }
        }

        template<typename Collision>
        void PostStep(Collision* collision, const site_t iFirstIndex, const site_t iSiteCount)
        {
          if (mVisControl->IsRendering())
          {
            collision->template DoPostStep<true>(iFirstIndex,
                                                 iSiteCount,
                                                 &mParams,
                                                 mLatDat,
                                                 propertyCache);
          }
          else
          {
            collision->template DoPostStep<false>(iFirstIndex,
                                                  iSiteCount,
                                                  &mParams,
                                                  mLatDat,
                                                  propertyCache);
          }
        }

        unsigned int inletCount;
        unsigned int outletCount;

        configuration::SimConfig *mSimConfig;
        net::Net* mNet;
        geometry::LatticeData* mLatDat;
        SimulationState* mState;
        iolets::BoundaryValues *mInletValues, *mOutletValues;

        LbmParameters mParams;
        vis::Control* mVisControl;

        const util::UnitConverter* mUnits;

        hemelb::reporting::Timers &timings;

        MacroscopicPropertyCache propertyCache;

        geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager;
    };

  } // Namespace lb
} // Namespace hemelb
#endif // HEMELB_LB_LB_H
