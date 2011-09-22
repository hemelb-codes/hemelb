// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

#include <math.h>
#include <limits>

#include "lb/lb.h"
#include "util/utilityFunctions.h"
#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace lb
  {
    hemelb::lb::LbmParameters *LBM::GetLbmParams()
    {
      return &mParams;
    }

    LBM::LBM(SimConfig *iSimulationConfig,
             net::Net* net,
             geometry::LatticeData* latDat,
             SimulationState* simState) :
      mSimConfig(iSimulationConfig), mNet(net), mLatDat(latDat), mState(simState),
          mParams(PULSATILE_PERIOD_s / (distribn_t) simState->GetTimeStepsPerCycle(),
                  latDat->GetVoxelSize())
    {
      ReadParameters();
    }

    void LBM::CalculateMouseFlowField(float densityIn,
                                      float stressIn,
                                      distribn_t &mouse_pressure,
                                      distribn_t &mouse_stress,
                                      double density_threshold_min,
                                      double density_threshold_minmax_inv,
                                      double stress_threshold_max_inv)
    {
      double density = density_threshold_min + densityIn / density_threshold_minmax_inv;
      double stress = stressIn / stress_threshold_max_inv;

      mouse_pressure = mUnits->ConvertPressureToPhysicalUnits(density * Cs2);
      mouse_stress = mUnits->ConvertStressToPhysicalUnits(stress);
    }

    void LBM::InitCollisions()
    {
      // TODO Note that the convergence checking is not yet implemented in the
      // new boundary condition hierarchy system.
      // It'd be nice to do this with something like
      // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

      kernels::InitParams initParams = kernels::InitParams();
      initParams.latDat = mLatDat;

      initParams.siteCount = mLatDat->GetInnerCollisionCount(0)
          + mLatDat->GetInterCollisionCount(0);
      mMidFluidCollision = new tMidFluidCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(1)
          + mLatDat->GetInterCollisionCount(1);
      mWallCollision = new tWallCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(2)
          + mLatDat->GetInterCollisionCount(2);
      initParams.boundaryObject = mInletValues;
      mInletCollision = new tInletOutletCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(3)
          + mLatDat->GetInterCollisionCount(3);
      initParams.boundaryObject = mOutletValues;
      mOutletCollision = new tInletOutletCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(4)
          + mLatDat->GetInterCollisionCount(4);
      initParams.boundaryObject = mInletValues;
      mInletWallCollision = new tInletOutletWallCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(5)
          + mLatDat->GetInterCollisionCount(5);
      initParams.boundaryObject = mOutletValues;
      mOutletWallCollision = new tInletOutletWallCollision(initParams);
    }

    void LBM::Initialise(site_t* iFTranslator,
                         vis::Control* iControl,
                         boundaries::BoundaryValues* iInletValues,
                         boundaries::BoundaryValues* iOutletValues,
                         util::UnitConverter* iUnits)
    {
      mInletValues = iInletValues;
      mOutletValues = iOutletValues;
      mUnits = iUnits;

      InitCollisions();

      receivedFTranslator = iFTranslator;

      SetInitialConditions();

      mVisControl = iControl;
    }

    void LBM::SetInitialConditions()
    {
      timeSpent = 0.0;

      distribn_t *f_old_p, *f_new_p, f_eq[D3Q15::NUMVECTORS];
      distribn_t density = 0.0;

      for (int i = 0; i < outlets; i++)
      {
        density += mOutletValues->GetDensityMin(i);
      }

      density /= outlets;

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        D3Q15::CalculateFeq(density, 0.0, 0.0, 0.0, f_eq);

        f_old_p = mLatDat->GetFOld(i * D3Q15::NUMVECTORS);
        f_new_p = mLatDat->GetFNew(i * D3Q15::NUMVECTORS);

        for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
        {
          f_new_p[l] = f_old_p[l] = f_eq[l];
        }
      }
    }

    void LBM::RequestComms()
    {
      timeSpent -= util::myClock();

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      for (std::vector<hemelb::topology::NeighbouringProcessor>::const_iterator it =
          netTop->NeighbouringProcs.begin(); it != netTop->NeighbouringProcs.end(); it++)
      {
        // Request the receive into the appropriate bit of FOld.
        mNet->RequestReceive<distribn_t> (mLatDat->GetFOld( (*it).FirstSharedF),
                                          (int) (*it).SharedFCount,
                                           (*it).Rank);

        // Request the send from the right bit of FNew.
        mNet->RequestSend<distribn_t> (mLatDat->GetFNew( (*it).FirstSharedF),
                                       (int) (*it).SharedFCount,
                                        (*it).Rank);

      }

      timeSpent += util::myClock();
    }

    void LBM::PreSend()
    {
      timeSpent -= util::myClock();

      site_t offset = mLatDat->GetInnerSiteCount();

      StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetInterCollisionCount(0));
      offset += mLatDat->GetInterCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, mLatDat->GetInterCollisionCount(1));
      offset += mLatDat->GetInterCollisionCount(1);

      mInletValues->FinishReceive();
      StreamAndCollide(mInletCollision, offset, mLatDat->GetInterCollisionCount(2));
      offset += mLatDat->GetInterCollisionCount(2);

      mOutletValues->FinishReceive();
      StreamAndCollide(mOutletCollision, offset, mLatDat->GetInterCollisionCount(3));
      offset += mLatDat->GetInterCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, mLatDat->GetInterCollisionCount(4));
      offset += mLatDat->GetInterCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetInterCollisionCount(5));

      timeSpent += util::myClock();
    }

    void LBM::PreReceive()
    {
      timeSpent -= util::myClock();

      site_t offset = 0;

      StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetInnerCollisionCount(0));
      offset += mLatDat->GetInnerCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, mLatDat->GetInnerCollisionCount(1));
      offset += mLatDat->GetInnerCollisionCount(1);

      StreamAndCollide(mInletCollision, offset, mLatDat->GetInnerCollisionCount(2));
      offset += mLatDat->GetInnerCollisionCount(2);

      StreamAndCollide(mOutletCollision, offset, mLatDat->GetInnerCollisionCount(3));
      offset += mLatDat->GetInnerCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, mLatDat->GetInnerCollisionCount(4));
      offset += mLatDat->GetInnerCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetInnerCollisionCount(5));

      timeSpent += util::myClock();
    }

    void LBM::PostReceive()
    {
      timeSpent -= util::myClock();

      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      for (site_t i = 0; i < netTop->TotalSharedFs; i++)
      {
        *mLatDat->GetFNew(receivedFTranslator[i])
            = *mLatDat->GetFOld(netTop->NeighbouringProcs[0].FirstSharedF + i);
      }

      // Do any cleanup steps necessary on boundary nodes
      site_t offset = 0;

      //TODO yup, this is horrible. If you read this, please improve the following code.
      PostStep(mMidFluidCollision, offset, mLatDat->GetInterCollisionCount(0));
      offset += mLatDat->GetInterCollisionCount(0);

      PostStep(mWallCollision, offset, mLatDat->GetInterCollisionCount(1));
      offset += mLatDat->GetInterCollisionCount(1);

      PostStep(mInletCollision, offset, mLatDat->GetInterCollisionCount(2));
      offset += mLatDat->GetInterCollisionCount(2);

      PostStep(mOutletCollision, offset, mLatDat->GetInterCollisionCount(3));
      offset += mLatDat->GetInterCollisionCount(3);

      PostStep(mInletWallCollision, offset, mLatDat->GetInterCollisionCount(4));
      offset += mLatDat->GetInterCollisionCount(4);

      PostStep(mOutletWallCollision, offset, mLatDat->GetInterCollisionCount(5));
      offset += mLatDat->GetInterCollisionCount(5);

      PostStep(mMidFluidCollision, offset, mLatDat->GetInnerCollisionCount(0));
      offset += mLatDat->GetInnerCollisionCount(0);

      PostStep(mWallCollision, offset, mLatDat->GetInnerCollisionCount(1));
      offset += mLatDat->GetInnerCollisionCount(1);

      PostStep(mInletCollision, offset, mLatDat->GetInnerCollisionCount(2));
      offset += mLatDat->GetInnerCollisionCount(2);

      PostStep(mOutletCollision, offset, mLatDat->GetInnerCollisionCount(3));
      offset += mLatDat->GetInnerCollisionCount(3);

      PostStep(mInletWallCollision, offset, mLatDat->GetInnerCollisionCount(4));
      offset += mLatDat->GetInnerCollisionCount(4);

      PostStep(mOutletWallCollision, offset, mLatDat->GetInnerCollisionCount(5));

      timeSpent += util::myClock();
    }

    void LBM::EndIteration()
    {
      timeSpent -= util::myClock();

      // Swap f_old and f_new ready for the next timestep.
      mLatDat->SwapOldAndNew();

      timeSpent += util::myClock();
    }

    // Update peak and average inlet velocities local to the current subdomain.
    void LBM::UpdateInletVelocities(unsigned long time_step)
    {
      distribn_t density;
      distribn_t vx, vy, vz;
      distribn_t velocity;

      int inlet_id;

      site_t offset = mLatDat->GetInnerCollisionCount(0) + mLatDat->GetInnerCollisionCount(1);

      for (site_t i = offset; i < offset + mLatDat->GetInnerCollisionCount(2); i++)
      {
        D3Q15::CalculateDensityAndVelocity(mLatDat->GetFOld(i * D3Q15::NUMVECTORS),
                                           density,
                                           vx,
                                           vy,
                                           vz);

        inlet_id = mLatDat->GetBoundaryId(i);

        vx *= inlet_normal[3 * inlet_id + 0];
        vy *= inlet_normal[3 * inlet_id + 1];
        vz *= inlet_normal[3 * inlet_id + 2];

        velocity = vx * vx + vy * vy + vz * vz;

        if (velocity > 0.)
        {
          velocity = sqrt(velocity) / density;
        }
        else
        {
          velocity = -sqrt(velocity) / density;
        }
      }

      offset = mLatDat->GetInnerSiteCount() + mLatDat->GetInterCollisionCount(0)
          + mLatDat->GetInterCollisionCount(1);

      for (site_t i = offset; i < offset + mLatDat->GetInterCollisionCount(2); i++)
      {
        D3Q15::CalculateDensityAndVelocity(mLatDat->GetFOld(i * D3Q15::NUMVECTORS),
                                           density,
                                           vx,
                                           vy,
                                           vz);

        inlet_id = mLatDat->GetBoundaryId(i);

        vx *= inlet_normal[3 * inlet_id + 0];
        vy *= inlet_normal[3 * inlet_id + 1];
        vz *= inlet_normal[3 * inlet_id + 2];

        velocity = vx * vx + vy * vy + vz * vz;

        if (velocity > 0.)
        {
          velocity = sqrt(velocity) / density;
        }
        else
        {
          velocity = -sqrt(velocity) / density;
        }
      }
    }

    /**
     * Return the amount of time spent doing lattice-Boltzmann
     * @return
     */
    double LBM::GetTimeSpent() const
    {
      return timeSpent;
    }

    // In the case of instability, this function restart the simulation
    // with twice as many time steps per period and update the parameters
    // that depends on this change.
    void LBM::Reset()
    {
      mState->DoubleTimeResolution();

      mParams.Update(PULSATILE_PERIOD_s / (distribn_t) mState->GetTimeStepsPerCycle(),
                     mLatDat->GetVoxelSize());

      SetInitialConditions();

      kernels::InitParams initParams = kernels::InitParams();
      initParams.latDat = mLatDat;

      initParams.siteCount = mLatDat->GetInnerCollisionCount(0)
          + mLatDat->GetInterCollisionCount(0);
      mMidFluidCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(1)
          + mLatDat->GetInterCollisionCount(1);
      mWallCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(2)
          + mLatDat->GetInterCollisionCount(2);
      initParams.boundaryObject = mInletValues;
      mInletCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(3)
          + mLatDat->GetInterCollisionCount(3);
      initParams.boundaryObject = mOutletValues;
      mOutletCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(4)
          + mLatDat->GetInterCollisionCount(4);
      initParams.boundaryObject = mInletValues;
      mInletWallCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(5)
          + mLatDat->GetInterCollisionCount(5);
      initParams.boundaryObject = mOutletValues;
      mOutletWallCollision->Reset(&initParams);
    }

    LBM::~LBM()
    {
      // Delete the translator between received location and location in f_new.
      delete[] receivedFTranslator;

      // Delete the collision and stream objects we've been using
      delete mMidFluidCollision;
      delete mWallCollision;
      delete mInletCollision;
      delete mOutletCollision;
      delete mInletWallCollision;
      delete mOutletWallCollision;

      // Delete various other arrays used
      delete[] inlet_normal;
    }
  }
}
