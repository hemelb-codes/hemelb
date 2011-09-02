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
    void LBM::RecalculateTauViscosityOmega()
    {
      mParams.Tau = 0.5 + (PULSATILE_PERIOD_s * BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3)
          / (Cs2 * ((double) mState->GetTimeStepsPerCycle() * mLatDat->GetVoxelSize()
              * mLatDat->GetVoxelSize()));

      mParams.Omega = -1.0 / mParams.Tau;
      mParams.StressParameter = (1.0 - 1.0 / (2.0 * mParams.Tau)) / sqrt(2.0);
      mParams.Beta = -1.0 / (2.0 * mParams.Tau); // This is just Omega / 2.0 by pure coincidence
    }

    hemelb::lb::LbmParameters *LBM::GetLbmParams()
    {
      return &mParams;
    }

    LBM::LBM(SimConfig *iSimulationConfig,
             net::Net* net,
             geometry::LatticeData* latDat,
             SimulationState* simState) :
      mSimConfig(iSimulationConfig), mNet(net), mLatDat(latDat), mState(simState)
    {
      // voxel_size = iSimulationConfig->VoxelSize;

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

      mouse_pressure = util::UnitConverter::ConvertPressureToPhysicalUnits(density * Cs2);
      mouse_stress = util::UnitConverter::ConvertStressToPhysicalUnits(stress);
    }

    template<typename tMidFluidCollision, typename tWallCollision, typename tInletOutletCollision,
        typename tInletOutletWallCollision, typename tCollisionOperator>
    void LBM::InitCollisions()
    {
      mStreamAndCollide
          = new hemelb::lb::collisions::StreamAndCollide<tMidFluidCollision, tWallCollision,
              tInletOutletCollision, tInletOutletWallCollision, tCollisionOperator>(mCollisionOperator);
      mPostStep = new hemelb::lb::collisions::PostStep<tMidFluidCollision, tWallCollision,
          tInletOutletCollision, tInletOutletWallCollision>();

      // TODO Note that the convergence checking is not yet implemented in the
      // new boundary condition hierarchy system.
      // It'd be nice to do this with something like
      // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

      mCollisions.resize(0);

      // WARNING: order is importnant
      mCollisions.push_back(new hemelb::lb::streamers::MidFluidCollision());
      mCollisions.push_back(new hemelb::lb::streamers::WallCollision());
      mCollisions.push_back(new hemelb::lb::streamers::InletOutletCollision(mInletValues));
      mCollisions.push_back(new hemelb::lb::streamers::InletOutletCollision(mOutletValues));
      mCollisions.push_back(new hemelb::lb::streamers::InletOutletWallCollision(mInletValues));
      mCollisions.push_back(new hemelb::lb::streamers::InletOutletWallCollision(mOutletValues));
    }

    void LBM::Initialise(site_t* iFTranslator,
                         vis::Control* iControl,
                         boundaries::BoundaryValues* iInletValues,
                         boundaries::BoundaryValues* iOutletValues)
    {
      mInletValues = iInletValues;

      mOutletValues = iOutletValues;

      mCollisionOperator = new CO(mLatDat, &mParams);

      InitCollisions<hemelb::lb::streamers::implementations::SimpleCollideAndStream<CO>,
          hemelb::lb::streamers::implementations::ZeroVelocityEquilibrium<CO>,
          hemelb::lb::streamers::implementations::NonZeroVelocityBoundaryDensity<CO>,
          hemelb::lb::streamers::implementations::ZeroVelocityBoundaryDensity<CO>, CO> ();

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
        density += mOutletValues->GetDensityMin(i);

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

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        // Wait for all boundaryComms to finish
        // TODO: This check is ugly, not ideal place or way to do it
        if (collision_type == 2)
        {
          mInletValues->FinishReceive();
        }
        else if (collision_type == 3)
        {
          mOutletValues->FinishReceive();
        }

        mCollisions[collision_type]->AcceptCollisionVisitor(mStreamAndCollide,
                                                            mVisControl->IsRendering(),
                                                            offset,
                                                            mLatDat->GetInterCollisionCount(collision_type),
                                                            &mParams,
                                                            mLatDat,
                                                            mVisControl);
        offset += mLatDat->GetInterCollisionCount(collision_type);
      }

      timeSpent += util::myClock();
    }

    void LBM::PreReceive()
    {
      timeSpent -= util::myClock();

      site_t offset = 0;

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        mCollisions[collision_type]->AcceptCollisionVisitor(mStreamAndCollide,
                                                            mVisControl->IsRendering(),
                                                            offset,
                                                            mLatDat->GetInnerCollisionCount(collision_type),
                                                            &mParams,
                                                            mLatDat,
                                                            mVisControl);
        offset += mLatDat->GetInnerCollisionCount(collision_type);
      }

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
      size_t offset = 0;

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        mCollisions[collision_type]->AcceptCollisionVisitor(mPostStep,
                                                            mVisControl->IsRendering(),
                                                            offset,
                                                            mLatDat->GetInnerCollisionCount(collision_type),
                                                            &mParams,
                                                            mLatDat,
                                                            mVisControl);
        offset += mLatDat->GetInnerCollisionCount(collision_type);
      }

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        mCollisions[collision_type]->AcceptCollisionVisitor(mPostStep,
                                                            mVisControl->IsRendering(),
                                                            offset,
                                                            mLatDat->GetInterCollisionCount(collision_type),
                                                            &mParams,
                                                            mLatDat,
                                                            mVisControl);
        offset += mLatDat->GetInterCollisionCount(collision_type);
      }

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

      RecalculateTauViscosityOmega();

      SetInitialConditions();

      mCollisionOperator->Reset(mLatDat, &mParams);
    }

    LBM::~LBM()
    {
      // Delete the translator between received location and location in f_new.
      delete[] receivedFTranslator;

      // Delete visitors
      delete mStreamAndCollide;
      delete mPostStep;

      // Delete Collision Operator
      delete mCollisionOperator;

      // Delete the collision and stream objects we've been using
      for (unsigned int i = 0; i < mCollisions.size(); i++)
      {
        delete mCollisions[i];
      }

      // Delete various other arrays used
      delete[] inlet_normal;
    }
  }
}
