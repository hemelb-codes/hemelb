#ifndef HEMELB_LB_H
#define HEMELB_LB_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "topology/NetworkTopology.h"
#include "lb/SimulationState.h"
#include "lb/kernels/Kernels.h"
#include "lb/collisions/Collisions.h"
#include "lb/streamers/Streamers.h"
#include "lb/boundaries/BoundaryValues.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"
#include "util/UnitConverter.h"
#include "SimConfig.h"
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
        //typedef kernels::LBGKNN<RHEO_MODEL> LB_KERNEL;

        // Standard LBGK collision operator
        typedef kernels::LBGK LB_KERNEL;

        typedef streamers::SimpleCollideAndStream<collisions::Normal<LB_KERNEL> >
            tMidFluidCollision;
        typedef streamers::SimpleCollideAndStream<collisions::ZeroVelocityEquilibrium<LB_KERNEL> >
            tWallCollision;
        typedef streamers::SimpleCollideAndStream<
            collisions::NonZeroVelocityEquilibriumFixedDensity<LB_KERNEL> > tInletOutletCollision;
        typedef streamers::SimpleCollideAndStream<collisions::ZeroVelocityEquilibriumFixedDensity<
            LB_KERNEL> > tInletOutletWallCollision;

      public:
        /**
         * Constructor, stage 1.
         * Object so initialized is not ready for simulation.
         * Must have Initialise(...) called also. Constructor separated due to need to access
         * the partially initialized LBM in order to initialize the arguments to the second construction phase.
         */
        LBM(hemelb::SimConfig *iSimulationConfig,
            net::Net* net,
            geometry::LatticeData* latDat,
            SimulationState* simState);
        ~LBM();


        void RequestComms(); ///< part of IteratedAction interface.
        void PreSend(); ///< part of IteratedAction interface.
        void PreReceive(); ///< part of IteratedAction interface.
        void PostReceive(); ///< part of IteratedAction interface.
        void EndIteration(); ///< part of IteratedAction interface.
        void Reset(); ///< part of IteratedAction interface.

        // TODO -- replace public member with accessor #26
        site_t total_fluid_sites;
        int inlets;

        // TODO -- replace built in type unsigned int with typedef #24
        void UpdateInletVelocities(unsigned long time_step); ///< Update peak and average inlet velocities local to the current subdomain.

        /**
         * Second constructor.
         *
         */
        void
        Initialise(site_t* iFTranslator,
                   vis::Control* iControl,
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
                                     PhysicalStress &mouse_stress
                                     );

        hemelb::lb::LbmParameters *GetLbmParams();
        double GetTimeSpent() const;

        site_t siteMins[3], siteMaxes[3];

      private:
        void SetInitialConditions();

        void InitCollisions();

        void ReadParameters();

        /***
         *  Calculate the BCs for each boundary site type and the
         *  non-equilibrium distribution functions.
         */
        void CalculateBC(distribn_t f[],
                         hemelb::geometry::LatticeData::SiteType const iSiteType,
                         unsigned int const iBoundaryId,
                         distribn_t *density,
                         distribn_t *vx,
                         distribn_t *vy,
                         distribn_t *vz,
                         distribn_t f_neq[]) const;

        void handleIOError(int iError);

        // Collision objects
        tMidFluidCollision* mMidFluidCollision;
        tWallCollision* mWallCollision;
        tInletOutletCollision* mInletCollision;
        tInletOutletCollision* mOutletCollision;
        tInletOutletWallCollision* mInletWallCollision;
        tInletOutletWallCollision* mOutletWallCollision;

        template<typename Collision>
        void StreamAndCollide(Collision* collision,
                              const site_t iFirstIndex,
                              const site_t iSiteCount)
        {
          if (mVisControl->IsRendering())
          {
            collision->template StreamAndCollide<true> (iFirstIndex,
                                                        iSiteCount,
                                                        &mParams,
                                                        mLatDat,
                                                        mVisControl);
          }
          else
          {
            collision->template StreamAndCollide<false> (iFirstIndex,
                                                         iSiteCount,
                                                         &mParams,
                                                         mLatDat,
                                                         mVisControl);
          }
        }

        template<typename Collision>
        void PostStep(Collision* collision, const site_t iFirstIndex, const site_t iSiteCount)
        {
          if (mVisControl->IsRendering())
          {
            collision->template DoPostStep<true> (iFirstIndex,
                                                  iSiteCount,
                                                  &mParams,
                                                  mLatDat,
                                                  mVisControl);
          }
          else
          {
            collision->template DoPostStep<false> (iFirstIndex,
                                                   iSiteCount,
                                                   &mParams,
                                                   mLatDat,
                                                   mVisControl);
          }
        }

        double timeSpent;

        double *inlet_normal;

        int outlets;

        SimConfig *mSimConfig;
        net::Net* mNet;
        geometry::LatticeData* mLatDat;
        SimulationState* mState;
        boundaries::BoundaryValues *mInletValues, *mOutletValues;

        LbmParameters mParams;
        vis::Control* mVisControl;

        util::UnitConverter* mUnits;

        site_t* receivedFTranslator;
    };
  }
}
#endif // HEMELB_LB_H
