#ifndef HEMELB_LB_H
#define HEMELB_LB_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "topology/NetworkTopology.h"
#include "lb/SimulationState.h"
#include "lb/collisions/Collisions.h"
#include "lb/collisions/CollisionVisitors.h"
#include "lb/collisions/implementations/Implementations.h"
#include "lb/collisions/implementations/CollisionOperators.h"
#include "lb/BoundaryComms.h"
#include "lb/BoundaryValues.h"
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
                   BoundaryComms* iInletComms,
                   BoundaryComms* iOutletComms,
                   util::UnitConverter* iUnits);

        void WriteConfigParallel(hemelb::lb::Stability stability,
                                 std::string output_file_name,
                                 BoundaryComms* iInletComms,
                                 BoundaryComms* iOutletComms);
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

        template<typename tMidFluidCollision, typename tWallCollision,
            typename tInletOutletCollision, typename tInletOutletWallCollision,
            typename tCollisionOperator>
        void InitCollisions(BoundaryComms* iInletComms, BoundaryComms* iOutletComms);

        void ReadParameters();

        void handleIOError(int iError);

        // Visitors
        hemelb::lb::collisions::CollisionVisitor* mStreamAndCollide;
        hemelb::lb::collisions::CollisionVisitor* mPostStep;

        // COllision Operator
        typedef hemelb::lb::collisions::implementations::LBGK CO;
        CO* mCollisionOperator;

        // Collision objects
        hemelb::lb::collisions::MidFluidCollision* mMidFluidCollision;
        hemelb::lb::collisions::WallCollision* mWallCollision;
        hemelb::lb::collisions::InletOutletCollision* mInletCollision;
        hemelb::lb::collisions::InletOutletCollision* mOutletCollision;
        hemelb::lb::collisions::InletOutletWallCollision* mInletWallCollision;
        hemelb::lb::collisions::InletOutletWallCollision* mOutletWallCollision;

        //TODO Get rid of this hack
        hemelb::lb::collisions::Collision* GetCollision(int i);

        double *inlet_normal;
        double voxel_size;
        int outlets;

        double mFileReadTime;

        SimConfig *mSimConfig;
        net::Net* mNet;
        geometry::LatticeData* mLatDat;
        SimulationState* mState;
        BoundaryValues* mBoundaryValues;

        LbmParameters mParams;
        vis::Control* mVisControl;

        util::UnitConverter* mUnits;

        site_t* receivedFTranslator;
    };
  }
}
#endif // HEMELB_LB_H
