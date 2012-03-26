#ifndef HEMELB_LB_LB_H
#define HEMELB_LB_LB_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "topology/NetworkTopology.h"
#include "lb/SimulationState.h"
#include "lb/kernels/Kernels.h"
#include "lb/kernels/momentBasis/MomentBases.h"
#include "lb/collisions/Collisions.h"
#include "lb/streamers/Streamers.h"
#include "lb/boundaries/BoundaryValues.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/UnitConverter.h"
#include "configuration/SimConfig.h"
#include "reporting/Timers.h"
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
    template<class LatticeType>
    class LBM : public net::IteratedAction
    {
      private:
        // TODO These should eventually be template parameters that are given to the generic LBM object.
        // At the moment, doing this will cause problems with other objects that have a pointer to the
        // LBM (and hence would also need to be templated, on the LBM-type).

        // Models of non-Newtonian rheology currently implemented.
        //typedef kernels::rheologyModels::CarreauYasudaRheologyModel RHEO_MODEL;
        //typedef kernels::rheologyModels::CassonRheologyModel RHEO_MODEL;
        //typedef kernels::rheologyModels::TruncatedPowerLawRheologyModel RHEO_MODEL;

        // LGBK operator with support for non-Newtonian flow
        //typedef kernels::LBGKNN<RHEO_MODEL, LatticeType> LB_KERNEL;

        // Multiple relaxation time collision operator
        //typedef kernels::MRT<kernels::momentBasis::DHumieresD3Q15MRTBasis> LB_KERNEL;

        // Standard LBGK collision operator
        typedef kernels::LBGK<LatticeType> LB_KERNEL;

        typedef streamers::SimpleCollideAndStream<collisions::Normal<LB_KERNEL> > tMidFluidCollision;
        typedef streamers::SimpleBounceBack<collisions::Normal<LB_KERNEL> > tWallCollision;
        typedef streamers::SimpleCollideAndStream<collisions::NonZeroVelocityEquilibriumFixedDensity<LB_KERNEL> > tInletOutletCollision;
        typedef streamers::SimpleCollideAndStream<collisions::ZeroVelocityEquilibriumFixedDensity<LB_KERNEL> > tInletOutletWallCollision;

      public:
        /**
         * Constructor, stage 1.
         * Object so initialized is not ready for simulation.
         * Must have Initialise(...) called also. Constructor separated due to need to access
         * the partially initialized LBM in order to initialize the arguments to the second construction phase.
         */
        LBM(hemelb::configuration::SimConfig *iSimulationConfig,
            net::Net* net,
            geometry::LatticeData* latDat,
            SimulationState* simState,
            reporting::Timer &atimer);
        ~LBM();

        void RequestComms(); ///< part of IteratedAction interface.
        void PreSend(); ///< part of IteratedAction interface.
        void PreReceive(); ///< part of IteratedAction interface.
        void PostReceive(); ///< part of IteratedAction interface.
        void EndIteration(); ///< part of IteratedAction interface.
        void Reset(); ///< part of IteratedAction interface.

        site_t TotalFluidSiteCount() const;
        void SetTotalFluidSiteCount(site_t);
        int InletCount() const {return inletCount;}
        int OutletCount() const {return outletCount;}

        /**
         * Second constructor.
         *
         */
        void Initialise(vis::Control* iControl,
                        boundaries::BoundaryValues* iInletValues,
                        boundaries::BoundaryValues* iOutletValues,
                        util::UnitConverter* iUnits);

        /**
         * This routine writes the flow field on file, using MPIO to coordinate
         * the writing. The format is detailed in io/formats/snapshot.h
         */
        // TODO filename argument should be const, but cannot be due to MPI constness issue #30
        void WriteConfigParallel(hemelb::lb::Stability const stability, std::string output_file_name) const;
        void ReadVisParameters();

        void CalculateMouseFlowField(const ScreenDensity densityIn,
                                     const ScreenStress stressIn,
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

        void ReadParameters();

        void handleIOError(int iError);

        // Collision objects
        tMidFluidCollision* mMidFluidCollision;
        tWallCollision* mWallCollision;
        tInletOutletCollision* mInletCollision;
        tInletOutletCollision* mOutletCollision;
        tInletOutletWallCollision* mInletWallCollision;
        tInletOutletWallCollision* mOutletWallCollision;

        template<typename Collision>
        void StreamAndCollide(Collision* collision, const site_t iFirstIndex, const site_t iSiteCount)
        {
          if (mVisControl->IsRendering())
          {
            collision->template StreamAndCollide<true>(iFirstIndex, iSiteCount, &mParams, mLatDat, propertyCache);
          }
          else
          {
            collision->template StreamAndCollide<false>(iFirstIndex, iSiteCount, &mParams, mLatDat, propertyCache);
          }
        }

        template<typename Collision>
        void PostStep(Collision* collision, const site_t iFirstIndex, const site_t iSiteCount)
        {
          if (mVisControl->IsRendering())
          {
            collision->template DoPostStep<true>(iFirstIndex, iSiteCount, &mParams, mLatDat, propertyCache);
          }
          else
          {
            collision->template DoPostStep<false>(iFirstIndex, iSiteCount, &mParams, mLatDat, propertyCache);
          }
        }

        unsigned int inletCount;
        unsigned int outletCount;

        configuration::SimConfig *mSimConfig;
        net::Net* mNet;
        geometry::LatticeData* mLatDat;
        SimulationState* mState;
        boundaries::BoundaryValues *mInletValues, *mOutletValues;

        LbmParameters mParams;
        vis::Control* mVisControl;

        util::UnitConverter* mUnits;

        reporting::Timer &timer;

        MacroscopicPropertyCache propertyCache;
    };


  } // Namespace lb
} // Namespace hemelb
#endif // HEMELB_LB_LB_H
