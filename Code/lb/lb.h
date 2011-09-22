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
#include "util/UnitConverter.h"
#include "vis/ColPixel.h"
#include "SimConfig.h"
#include <typeinfo>

namespace hemelb
{
  namespace lb
  {
    class LBM : public net::IteratedAction
    {
      private:
        // TODO These should eventually be template parameters that are given to the generic LBM object.
        // At the moment, doing this will cause problems with other objects that have a pointer to the
        // LBM (and hence would also need to be templated, on the LBM-type).
        typedef streamers::SimpleCollideAndStream<collisions::Normal<kernels::LBGK> >
            tMidFluidCollision;
        typedef streamers::SimpleCollideAndStream<
            collisions::ZeroVelocityEquilibrium<kernels::LBGK> > tWallCollision;
        typedef streamers::SimpleCollideAndStream<
            collisions::NonZeroVelocityEquilibriumFixedDensity<kernels::LBGK> >
            tInletOutletCollision;
        typedef streamers::SimpleCollideAndStream<collisions::ZeroVelocityEquilibriumFixedDensity<
            kernels::LBGK> > tInletOutletWallCollision;

      public:
        LBM(hemelb::SimConfig *iSimulationConfig,
            net::Net* net,
            geometry::LatticeData* latDat,
            SimulationState* simState);
        ~LBM();

        void RequestComms();
        void PreSend();
        void PreReceive();
        void PostReceive();
        void EndIteration();
        void Reset();

        site_t total_fluid_sites;
        int inlets;

        void UpdateBoundaryDensities(unsigned long time_step);
        void UpdateInletVelocities(unsigned long time_step);

        void
        Initialise(site_t* iFTranslator,
                   vis::Control* iControl,
                   boundaries::BoundaryValues* iInletValues,
                   boundaries::BoundaryValues* iOutletValues,
                   util::UnitConverter* iUnits);

        void WriteConfigParallel(hemelb::lb::Stability stability, std::string output_file_name);
        void ReadVisParameters();

        void CalculateMouseFlowField(float densityIn,
                                     float stressIn,
                                     distribn_t &mouse_pressure,
                                     distribn_t &mouse_stress,
                                     double density_threshold_min,
                                     double density_threshold_minmax_inv,
                                     double stress_threshold_max_inv);

        hemelb::lb::LbmParameters *GetLbmParams();
        double GetTimeSpent() const;

        site_t siteMins[3], siteMaxes[3];

      private:
        void SetInitialConditions();

        void InitCollisions();

        void ReadParameters();
        void CalculateBC(distribn_t f[],
                         hemelb::geometry::LatticeData::SiteType iSiteType,
                         unsigned int iBoundaryId,
                         distribn_t *density,
                         distribn_t *vx,
                         distribn_t *vy,
                         distribn_t *vz,
                         distribn_t f_neq[]);

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
            collision->template StreamAndCollide<true>(iFirstIndex,
                                                       iSiteCount,
                                                       &mParams,
                                                       mLatDat,
                                                       mVisControl);
          }
          else
          {
            collision->template StreamAndCollide<false>(iFirstIndex,
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
