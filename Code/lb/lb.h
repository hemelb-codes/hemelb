#ifndef HEMELB_LB_H
#define HEMELB_LB_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "topology/NetworkTopology.h"
#include "lb/SimulationState.h"
#include "lb/streamers/Collisions.h"
#include "lb/collisions/CollisionVisitors.h"
#include "lb/streamers/Implementations.h"
#include "lb/collisions/CollisionOperators.h"
#include "vis/ColPixel.h"
#include "SimConfig.h"
#include <typeinfo>

namespace hemelb
{
  namespace lb
  {
    class LBM : public net::IteratedAction
    {
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

        distribn_t ConvertPressureToLatticeUnits(double pressure) const;
        distribn_t ConvertVelocityToLatticeUnits(double velocity) const;
        distribn_t ConvertStressToLatticeUnits(double stress) const;

        void UpdateBoundaryDensities(unsigned long time_step);
        void UpdateInletVelocities(unsigned long time_step);

        void Initialise(site_t* iFTranslator, vis::Control* iControl);

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

        site_t siteMins[3], siteMaxes[3];

      private:
        void RecalculateTauViscosityOmega();
        void SetInitialConditions();

        double ConvertPressureToPhysicalUnits(double distribn_t) const;
        double ConvertStressToPhysicalUnits(double distribn_t) const;
        double ConvertVelocityToPhysicalUnits(double distribn_t) const;

        void CalculateBC(distribn_t f[],
                         hemelb::geometry::LatticeData::SiteType iSiteType,
                         unsigned int iBoundaryId,
                         distribn_t *density,
                         distribn_t *vx,
                         distribn_t *vy,
                         distribn_t *vz,
                         distribn_t f_neq[]);

        template<typename tMidFluidCollision, typename tWallCollision,
            typename tInletOutletCollision, typename tInletOutletWallCollision>
        void InitCollisions();

        void ReadParameters();

        void allocateInlets(int nInlets);
        void allocateOutlets(int nOutlets);

        void handleIOError(int iError);

        distribn_t ConvertPressureGradToLatticeUnits(double pressure_grad) const;
        double ConvertPressureGradToPhysicalUnits(distribn_t pressure_grad) const;

        // Visitors
        hemelb::lb::collisions::CollisionVisitor* mStreamAndCollide;
        hemelb::lb::collisions::CollisionVisitor* mPostStep;

        // Collision objects
        hemelb::lb::streamers::MidFluidCollision* mMidFluidCollision;
        hemelb::lb::streamers::WallCollision* mWallCollision;
        hemelb::lb::streamers::InletOutletCollision* mInletCollision;
        hemelb::lb::streamers::InletOutletCollision* mOutletCollision;
        hemelb::lb::streamers::InletOutletWallCollision* mInletWallCollision;
        hemelb::lb::streamers::InletOutletWallCollision* mOutletWallCollision;

        //TODO Get rid of this hack
        hemelb::lb::streamers::Collision* GetCollision(int i);

        distribn_t *inlet_density_avg, *inlet_density_amp;
        distribn_t *outlet_density_avg, *outlet_density_amp;
        distribn_t *inlet_density_phs, *outlet_density_phs;
        distribn_t* inlet_density, *outlet_density;
        double *inlet_normal;
        double voxel_size;
        int outlets;

        double mFileReadTime;

        SimConfig *mSimConfig;
        net::Net* mNet;
        geometry::LatticeData* mLatDat;
        SimulationState* mState;

        LbmParameters mParams;
        vis::Control* mVisControl;

        site_t* receivedFTranslator;
    };
  }
}
#endif // HEMELB_LB_H
